import torch
import numpy as np
import time
from torch.utils.data import random_split


def mask_spectra_data(spectra, spectra_mask, remain_ratio=0.1, spectra_zero_ratio=0.1, device='cpu'):
    """Args:
            spectra: float Tensor (batch, 300, 2)
            spectra_mask: float Tensor (batch, 300)
            remain_ratio: float=0.1 (ratio of samples that are not masked)
            spectra_zero_ratio: float=0.1 (per-sample spectrum masking ratio)
            device: cpu or cuda
        Returns:
            spectra: float Tensor (batch, 300, 2)
            spectra_mask: float Tensor (batch, 300)
    """
    # Randomly remove 10% of <mz, intensity> pairs
    spectra_mask_backup = torch.rand(spectra.shape[0], spectra.shape[1]).to(device)
    spectra_mask_backup = torch.where(spectra_mask_backup < spectra_zero_ratio, 0, 1)

    # Keep 10% of spectra unchanged
    full_indices = list(range(spectra.shape[0]))
    remain_size = int(remain_ratio * len(full_indices))
    mask_size = len(spectra) - remain_size
    remain_spectra_index_list, mask_spectra_index_list = random_split(full_indices, [remain_size, mask_size])

    # Set mask values to 1 for the "remain" samples (unchanged)
    for i in remain_spectra_index_list:
        spectra_mask_backup[i] = 1

    # Apply the mask matrix to spectrum and spectra_mask
    spectra = torch.mul(spectra, torch.unsqueeze(spectra_mask_backup, dim=-1)).to(device)
    spectra_mask = torch.mul(spectra_mask, spectra_mask_backup).to(device)
    return spectra, spectra_mask


def build_contiguous_chunks_mask_simple(
    lengths: torch.Tensor,     # (B_sel,)
    total_mask_ratio: float,   # e.g. 0.3
    max_chunks: int = 3,
    device: str | torch.device = "cpu",
) -> torch.Tensor:
    """
    Generate 1..max_chunks non-overlapping contiguous chunk masks for each sample.
    - For each sample, randomly choose n_chunks ∈ {1..max_chunks}
    - Each chunk length ≈ L * total_mask_ratio / n_chunks (floored, at least 1)
    - Randomly choose a start position from currently unmasked positions
    Returns: (B_sel, L_max) bool
    """
    device = torch.device(device)
    lengths = lengths.to(device)                # (B_sel,)
    B_sel = lengths.shape[0]
    if B_sel == 0:
        return torch.zeros(0, 0, dtype=torch.bool, device=device)

    # Only generate masks for samples with length > 0; rows with length 0 remain all-False
    valid = lengths > 0                         # (B_sel,)
    if not valid.any():
        return torch.zeros(B_sel, 0, dtype=torch.bool, device=device)

    L_max = int(lengths[valid].max().item())
    if L_max <= 0:
        return torch.zeros(B_sel, 0, dtype=torch.bool, device=device)

    # Random number of chunks (1..max_chunks)
    n_chunks = torch.randint(1, max_chunks + 1, (B_sel,), device=device)   # (B_sel,)

    # Base chunk length per sample (compute in float first, then cast to long)
    lengths_f = lengths.float()
    base_chunk_len = (lengths_f * total_mask_ratio / n_chunks.float()).floor().long()  # (B_sel,)
    # Ensure >= 1
    base_chunk_len = torch.maximum(base_chunk_len, torch.ones_like(base_chunk_len))
    # Ensure <= length
    base_chunk_len = torch.minimum(base_chunk_len, lengths)                             # (B_sel,)

    mask = torch.zeros(B_sel, L_max, dtype=torch.bool, device=device)

    # A small loop over chunk index; each iteration is fully vectorized over the batch dimension
    for k in range(max_chunks):
        # Only for samples that need the k-th chunk and have length > 0
        active = (n_chunks > k) & valid                      # (B_sel,)
        if not active.any():
            continue

        L_active = lengths[active]                           # (B_act,)
        chunk_len_active = base_chunk_len[active]            # (B_act,)

        current_mask = mask[active]                          # (B_act, L_max)
        pos = torch.arange(L_max, device=device).view(1, -1) # (1, L_max)
        within_len = pos < L_active.unsqueeze(1)             # (B_act, L_max)
        free_pos = within_len & (~current_mask)              # (B_act, L_max)

        # If the length is smaller than chunk_len, we cannot place the chunk; skip those samples
        can_place = L_active >= chunk_len_active             # (B_act,)
        if not can_place.any():
            continue

        L_act2 = L_active[can_place]                         # (B_act2,)
        len_act2 = chunk_len_active[can_place]               # (B_act2,)

        # max_start = max(L_act2 - len_act2, 0)
        max_start = L_act2 - len_act2                        # (B_act2,)
        max_start = torch.maximum(max_start, torch.zeros_like(max_start))

        rand_u = torch.rand_like(max_start.float())
        start = (rand_u * (max_start + 1).float()).floor().long()  # (B_act2,)
        end = start + len_act2                                     # (B_act2,)
        # Defensive clamp: end <= L_act2
        end = torch.minimum(end, L_act2)

        pos2 = pos                                                 # (1, L_max)
        in_chunk2 = (pos2 >= start.unsqueeze(1)) & (pos2 < end.unsqueeze(1))  # (B_act2, L_max)

        free_pos2 = free_pos[can_place]
        in_chunk2 = in_chunk2 & free_pos2

        # Write back into mask[active & can_place]
        submask = mask[active]
        submask[can_place] |= in_chunk2
        mask[active] = submask

    return mask


def mask_batch_with_mass_mlm(
    batch,
    pad_index: int = 0,
    mask_index: int = 1,
    unk_index: int = 28,              # keep the signature; no longer used internally
    mlm_total_mask_ratio: float = 0.3,
    mlm_sample_ratio: float = 1.0,    # do MLM for the whole batch by default
    mlm_max_chunks: int = 3,
    spectrum_mask_ratio: float = 0.9,
    spectrum_drop_ratio: float = 0.1,
    device: str = "cpu",
):
    """
    Updated version:

    - Apply contiguous token masking over the full batch
    - The same batch participates in both DDA and MLM training
    - Masked positions are replaced with the <mask> token (mask_index); unmasked positions keep the original token

    Returns:
        spectra, spectra_mask, precursors,
        tokens_input, tokens_label_out, token_mask_out,
        label, weight
    """
    spectra, spectra_mask, precursors, tokens, label, weight = batch
    device = torch.device(device)

    spectra = spectra.to(device).float()
    spectra_mask = spectra_mask.to(device).bool()
    precursors = precursors.to(device).float()
    tokens = tokens.to(device).long()
    label = label.to(device).float()
    weight = weight.to(device).float()

    B, L = tokens.shape

    # ========== 1. Spectrum augmentation ==========
    spectra, spectra_mask = mask_spectra_data(
        spectra,
        spectra_mask,
        remain_ratio=1 - spectrum_mask_ratio,
        spectra_zero_ratio=spectrum_drop_ratio,
        device=device,
    )

    # ========== 2. MLM masking over the full batch ==========
    rand_for_mlm = torch.rand(B, device=device)
    mlm_mask_batch = rand_for_mlm < mlm_sample_ratio      # (B,)

    valid_mask = (tokens != pad_index)                    # (B, L)
    lengths = valid_mask.sum(dim=1)                       # (B,)

    tokens_input = tokens.clone()
    tokens_label_out = torch.zeros_like(tokens, dtype=torch.long, device=device)
    token_mask_out = torch.zeros_like(tokens, dtype=torch.bool, device=device)

    # Only generate chunk masks for samples that do MLM and have length > 0
    idx_M_all = torch.nonzero(mlm_mask_batch, as_tuple=False).view(-1)     # (B_M,)
    has_len = lengths > 0
    if idx_M_all.numel() > 0:
        idx_M_valid = idx_M_all[has_len[idx_M_all]]
        if idx_M_valid.numel() > 0:
            lengths_M = lengths[idx_M_valid]  # (B_M_valid,)

            # Generate contiguous chunk masks (tokens only)
            chunk_mask_M = build_contiguous_chunks_mask_simple(
                lengths_M,
                total_mask_ratio=mlm_total_mask_ratio,
                max_chunks=mlm_max_chunks,
                device=device,
            )  # (B_M_valid, L_max)

            L_max = chunk_mask_M.shape[1]
            if L_max < L:
                pad_cols = L - L_max
                pad = torch.zeros(chunk_mask_M.shape[0], pad_cols, dtype=torch.bool, device=device)
                chunk_mask_M = torch.cat([chunk_mask_M, pad], dim=1)
            elif L_max > L:
                chunk_mask_M = chunk_mask_M[:, :L]

            # Only valid tokens (non-pad) can be masked
            valid_mask_M = valid_mask[idx_M_valid]          # (B_M_valid, L)
            chunk_mask_M = chunk_mask_M & valid_mask_M      # (B_M_valid, L)

            # === Key change: pure MLM on masked positions, no mass token ===
            # tokens_label_out: store original tokens only at masked positions (MLM targets)
            tokens_M = tokens[idx_M_valid]                  # (B_M_valid, L)
            tokens_label_out[idx_M_valid] = torch.where(
                chunk_mask_M,
                tokens_M,
                torch.zeros_like(tokens_M),
            ).long()

            # tokens_input: masked positions -> <mask>, others keep original tokens
            tokens_input[idx_M_valid] = torch.where(
                chunk_mask_M,
                torch.full_like(tokens_M, mask_index, dtype=torch.long),
                tokens_M,
            )

            # token_mask_out: mark MLM positions only (masked positions)
            token_mask_out[idx_M_valid] = chunk_mask_M

    return (
        spectra,
        spectra_mask,
        precursors,
        tokens_input,
        tokens_label_out,
        token_mask_out,
        label,
        weight,
    )


def mask_batch_with_augmentation(
    batch,
    pad_index: int = 0,
    mask_index: int = 1,
    unk_index: int = 28,
    token_mask_ratio: float = 0.3,
    mlm_sample_ratio: float = 0.9,
    mlm_max_chunks: int = 3,
    spectrum_topk: int = 300,
    spectrum_frac_shuffle: float = 0.7,
    spectrum_num_shuffle: int = 3,
    spectrum_mask_ratio: float = 0.9,
    spectrum_drop_ratio: float = 0.1,
    spectrum_augment: bool = True,
    device: str = "cpu",
):
    (
        spectra,
        spectra_mask,
        precursors,
        tokens_input,
        tokens_label,
        token_mask,
        label,
        weight,
    ) = mask_batch_with_mass_mlm(
        batch,
        pad_index=pad_index,
        mask_index=mask_index,
        unk_index=unk_index,
        mlm_total_mask_ratio=token_mask_ratio,
        mlm_sample_ratio=mlm_sample_ratio,
        mlm_max_chunks=mlm_max_chunks,
        spectrum_mask_ratio=spectrum_mask_ratio,
        spectrum_drop_ratio=spectrum_drop_ratio,
        device=device,
    )

    return (
        spectra,
        spectra_mask,
        precursors,
        tokens_input,
        tokens_label,
        token_mask,
        label,
        weight,
    )
