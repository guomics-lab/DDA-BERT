<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://github.com/zhiyuajun/DDA-BERT/blob/main/DDA-BERT.png" width="240" height="240">
</p>

# Description
DDA-BERT is an end-to-end deep learning tool for rescoring peptide-spectrum matches (PSMs) in data-dependent acquisition (DDA) proteomics. Built on a Transformer-based architecture and trained on 12,285 DDA-MS files comprising approximately 271 million high-confidence PSMs, DDA-BERT effectively models the complex relationships between peptide sequences and tandem mass spectra. It supports rescoring results generated from a single database search engine as well as integrated outputs from multiple search engines.
DDA-BERT demonstrates robust and consistent performance across diverse biological systems, including animal, plant, and microbial proteomes. When applied to HLA immunopeptidomics datasets, DDA-BERT continues to show strong and reliable performance. It is applicable to low-input contexts such as trace-level and single-cell proteomics, offering a scalable and reliable solution for improving peptide identification in mass spectrometry-based workflows.

# Installation
DDA-BERT runs without installation and does not require any additional environment configuration. The executable is available at https://guomics.com/software/DDA-BERT.

Hardware Requirements:  
•	Operating System: Compatible with Linux-based operating systems.  
•	Processor: A dual-core processor is recommended; the platform can also run on a single-core processor.  
•	Memory: At least 40 GB of RAM is recommended. Higher memory configurations are advised for processing large-scale mass spectrometry or FASTA datasets.  
•	Storage: A minimum of 100 GB of available disk space is recommended.  
•	Graphics Processing Unit (GPU): An NVIDIA GPU that supports bfloat16 (bf16) precision inference is required. CUDA support is necessary, and a minimum of 20GB GPU memory is recommended.

## Directories Overview
This repository is organized into several key directories:

### Software
The `software` directory contains the source code for the identification software along with instructions on how to use it. Here, you will find all necessary scripts and binaries required for running the software.

### Training and Evaluation
Under the `training_eval` directory, you will discover the code used for training and evaluating models. This includes scripts for preparing datasets, configuring training parameters, and executing training jobs.

Scripts
The `scripts` directory holds various scripts aimed at facilitating the process of plotting and visualizing data. These tools can be particularly useful for generating insights from experimental results.

### Demo Data
For demonstration purposes, we have included a set of sample data in the `demo_data` directory. This data can be used to test the functionalities of the software and models without needing to prepare your own dataset.


## Evaluation

Evaluation typically completes in about 20 minutes, depending on the number of spectra and available GPU/CPU resources.

**For the provided test data, the complete workflow, including database search, PSM rescoring, and protein inference, takes approximately 17.5 minutes on a single NVIDIA A100 (40GB) GPU with 20 CPU cores (AMD EPYC 7742 64-Core Processor).**

Currently, the tool supports only the .mzML format. However, the full source code is openly available and modifiable (see license for details), allowing users to adapt the tool to accommodate other data formats as needed. Future versions will gradually introduce direct compatibility with additional commonly used mass spectrometry formats, such as Sciex .wiff, Bruker .d, and other raw data types.

## Results
Benchmarking results: https://zenodo.org/records/15923904

Results are output in CSV format as a comprehensive summary table that is easy to manipulate and interpret, facilitating further biological insights and downstream applications.

# License
This software is licensed under a custom license that allows personal use but prohibits commercial use. For more details, see the LICENSE file.

# Contact
For any questions or licensing inquiries, please contact: Dr. Guo E-mail: guotiannan@westlake.edu.cn
www.guomics.com
