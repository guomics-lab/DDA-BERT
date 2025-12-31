## Description
DDA-BERT: a computational platform designed for data-dependent acquisition (DDA) proteomics analysis.

## Software
The DDA-BERT executable and user manual are available at https://guomics.com/software/DDA-BERT. Run it directly without installation or additional configuration.

# Quick Start Guide

## Installation

DDA-BERT can be installed and used on Linux systems.

There are different types of installation or use possible:

- [**Portable Executable:**](#portable-executable) Choose this ready-to-run, no-install version!
- [**Developer installation:**](#developer-installation) Choose this installation if you
  are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/)
  and Python. This installation allows access to all available features
  of DDA-BERT and even allows to modify its source code directly. You can use DDA-BERT as a Python package with this choose.
- [**Docker installation**](#docker-installation) Choose this installation if you want to use DDA-BERT with docker.

### Portable Executable

You can download the latest release of DDA-BERT [here](https://guomics.com/software/DDA-BERT/downloads.html).

### Developer installation

#### Install with pip
It is strongly recommended to install DDA-BERT in its own environment. 
1. Open the console and create a new conda environment: conda create --name dda-bert python=3.10 
2. Activate the environment: conda activate dda-bert 
3. Install DDA-BERT via pip: pip install dda-bert.

#### Install with Developer from source
⏱️ Estimated setup time: **~2–5 minutes**

#### Prerequisites
Ensure that Conda or Miniconda is installed. We recommend Miniconda; follow the installation instructions on the official website.

```shell
git clone https://github.com/guomics-lab/DDA-BERT.git
cd software/DDA-BERT
```

```shell
conda create -n DDA-BERT python=3.10
conda activate DDA-BERT
```

```shell
pip install uv
uv pip install -e . --refresh
```

Linux command-line run
```shell
dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output-path=/out/
dda-bert score --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --sage-file-dir=xxx --fp-file-dir=xxx --ap-file-dir=xxx --engines=sage,fp,ap  --output-path=/out/
```


### Docker installation


This repository provides a self-contained **Docker image** that encapsulates all necessary environments and
dependencies for the MassNet-DDA conversion utility. By using this image, users can quickly launch the tool without
complex setup.

### Prerequisites

* **Docker Engine** (for Linux) must be installed and running.

---

The process involves **pulling the image from Docker Hub** and then running a container, mapping your local data
directory to the container's working directory.

1. **Ensure the Docker service is running.**
2. **Pull the Docker Image** from the registry in your terminal:
   ```bash
   docker pull guomics2017/dda-bert:v3.0
   ```
3. **Run the Container** (Example using a typical Linux absolute path):
   ```bash
   docker run --rm -v /home/user/DDA-BERT:/home/test_data guomics2017/dda-bert:v3.0 assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output-path=/out/
   ```


## Hardware Requirements:
•	Operating System: Compatible with Linux-based operating systems.

•	Processor: A dual-core processor is recommended; the platform can also run on a single-core processor.

•	Memory: At least 40 GB of RAM is recommended. Higher memory configurations are advised for processing large-scale mass spectrometry or FASTA datasets.

•	Graphics Processing Unit (GPU): An NVIDIA GPU that supports bfloat16 (bf16) precision inference is required. CUDA support is necessary, and a minimum of 20GB GPU memory is recommended.

## License
This software is licensed under a custom license that allows academic use but prohibits commercial use. For more details, see the LICENSE file.

## Contact
For any questions or licensing inquiries, please contact:
Dr Guo
E-mail: guotiannan@westlake.edu.cn
www.guomics.com