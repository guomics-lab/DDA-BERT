# DDA-BERT
# Description
DDA-BERT is an open-source, end-to-end deep learning tool for rescoring peptide-spectrum matches (PSMs) in data-dependent acquisition (DDA) proteomics. Built on a Transformer-based architecture and trained on 3,701 DDA-MS files encompassing approximately 82 million high-confidence PSMs, it effectively models the complex relationships between peptide sequences and tandem mass spectra. DDA-BERT demonstrates robust and consistent performance across a diverse range of biological systems, including animal, plant, and microbial proteomes. It is particularly effective in low-input contexts such as trace-level and single-cell proteomics, offering a scalable and reliable solution for improving peptide identification in mass spectrometry-based workflows.

Currently, the tool supports only the .mzML format. However, the full source code is openly available and modifiable (see license for details), allowing users to adapt the tool to accommodate other data formats as needed. Future versions will gradually introduce direct compatibility with additional commonly used mass spectrometry formats, such as Sciex .wiff, Bruker .d, and other raw data types.

# Installation
DDA-BERT runs install-free and requires no additional configuration of the environment. 

Hardware Requirements:  
•	Operating System: Compatible with Linux-based operating systems.  
•	Processor: A dual-core processor is recommended; the platform can also run on a single-core processor.  
•	Memory: At least 40 GB of RAM is recommended. Higher memory configurations are advised for processing large-scale mass spectrometry or FASTA datasets.  
•	Storage: A minimum of 100 GB of available disk space is recommended.  
•	Graphics Processing Unit (GPU): An NVIDIA GPU that supports bfloat16 (bf16) precision inference is required. CUDA support is necessary, and a minimum of 20GB GPU memory is recommended.

The DDA-BERT executable is available at https://guomics.com/software/DDA-BERT and on Zenodo (https://zenodo.org/records/15923904).

## Run Instructions
## Step1: Download Model and Test Files
You may also use your own .raw and .mzML files.

Model checkpoint: ./software/resource/model/mp_rank_00_model_states.pt 

Alternatively, it can be downloaded from Zenodo: https://zenodo.org/records/15923904

Test data files: ./demo_data/HeLa_digest_SPME_1ng_1.mzML and ./demo_data/HeLa_digest_SPME_1ng_1.raw

##Note: .raw files can be converted to .mzML format using the MSConvertGUI tool from ProteoWizard. The default settings are sufficient, or you may refer to the configuration file DB_search_config/msConvert.config.txt for custom conversion options.
During execution, make sure that the .mzML and corresponding .raw file are placed in the same directory.

## Step2: Run the Command
Execute the following command in your terminal:
```shell
cd DDA-BERT; 
./DDA-BERT --mzml_paths=/data/example.mzML --fasta=/data/example.fasta --output_path=/out/
```

##The evaluation typically takes approximately 20 minutes, depending on the number of spectra, as well as the number of GPUs and CPU cores available. For reference, using a single NVIDIA A100 (40GB) GPU and 20 CPU cores, the complete workflow—including database search, PSM rescoring, and protein inference—was completed in about 23.5 minutes.


##The evaluation typically takes approximately 20 minutes, depending on the number of spectra, as well as the number of GPUs and CPU cores available. For reference, using a single NVIDIA A100 (40GB) GPU and 20 CPU cores, the complete workflow—including database search, PSM rescoring, and protein inference—was completed in about 23.5 minutes.

## Results
Benchmarking results: https://zenodo.org/records/15923904

Results are output in CSV format as a comprehensive summary table that is easy to manipulate and interpret, facilitating further biological insights and downstream applications.

# License
This software is licensed under a custom license that allows personal use but prohibits commercial use. For more details, see the LICENSE file.

# Contact
For any questions or licensing inquiries, please contact: Dr Guo E-mail: guotiannan@westlake.edu.cn
www.guomics.com
