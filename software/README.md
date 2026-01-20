# Quick Start Guide

## Installation

DDA-BERT is designed to run on Linux systems and can be deployed using multiple installation options to accommodate different usage scenarios.

Users may choose one of the following installation or usage modes:

- [**Portable Executable:**](#portable-executable) A ready-to-run, no-install option recommended for quick evaluation and standard usage.
- [**Python Package Installation:**](#python-package-installation) This option installs DDA-BERT via pip and enables users to modify the source code as needed, for example, to support rescoring outputs from additional database search engines. In this mode, DDA-BERT can be used directly as a Python package within custom workflows.
- [**Docker installation**](#docker-installation) Recommended for users who prefer containerized environments or require reproducible deployments.

### Portable Executable

#### Step1: Download

Download the latest DDA-BERT portable executable and accompanying test files from the official release page.

You can download the latest release of DDA-BERT [https://guomics.com/software/DDA-BERT](https://guomics.com/software/DDA-BERT/downloads.html).

#### Step2: Run
Unzip the downloaded archive and execute the following command in a terminal:
```shell
cd DDA-BERT; 
./dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output-path=/out/
```

### Installation from Source (Editable Mode)

#### Python Package Installation
It is strongly recommended to install DDA-BERT in an isolated Conda environment.
⏱️ Estimated setup time: **~10–15 minutes**

**Test data:** demo_data/HeLa_digest_SPME_1ng_1.mzML and demo_data/HeLa_digest_SPME_1ng_1.raw. You can also use your own .raw and .mzML files.

#### Prerequisites

Download jdk11 from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/DDA-BERT/jdk-11.0.26.zip),
unzip and move to project root directory.
> **⚠️Note**: This project was built using **FragPipe v22**, which includes the following core components: **MSFragger v4.1**, **Philosopher v5.1.1**, **diaTracer v1.1.5**, **IonQuant v1.10.27**. If you're using a different version of FragPipe, make sure to download compatible versions of these tools and properly configure your environment to ensure smooth execution of the analysis pipeline.

You can download FragPipe from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/DDA-BERT/FragPipe22_0.zip) with the following core components.


#### Prerequisites
Clone the DDA-BERT repository and set up an isolated Conda environment:

```shell
git clone https://github.com/guomics-lab/DDA-BERT.git
cd software
```

```shell
conda create -n DDA-BERT python=3.10
conda activate DDA-BERT
```

```shell
pip install uv
uv pip install -e . --refresh
```

Running DDA-BERT (Linux Command Line)

PSM Assessment
```shell
dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output-path=/out/
```

Multi-engine PSM Rescoring
```shell
dda-bert score --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --sage-file-dir=xxx --fp-file-dir=xxx --ap-file-dir=xxx --engines=sage,fp,ap  --output-path=/out/
```
### Docker Installation

DDA-BERT is available as a self-contained Docker image that includes all required dependencies and runtime environments. This option enables users to run DDA-BERT without manual environment configuration and is well suited for reproducible and portable deployments.

### Prerequisites

* **Docker Engine** installed and running on a Linux system
* Sufficient permissions to pull images from Docker Hub and run containers

Ensure that the Docker service is active before proceeding.

Step 1: Pull the Docker Image
Pull the pre-built DDA-BERT image from Docker Hub:
   ```bash
   docker pull guomics2017/dda-bert:v3.1
   ```
Step 2: Run the Container

Below is an example command using a typical Linux absolute path, where a local directory is mounted into the container for data access:

   ```bash
   docker run --rm -v /home/user/DDA-BERT:/home/test_data guomics2017/dda-bert:v3.1 assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output-path=/out/
   ```
In this example, the local directory /home/user/DDA-BERT is mounted into the container and used as the working directory for input and output files.
