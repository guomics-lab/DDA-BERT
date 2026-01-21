# Quick Start Guide

Multiple installation options are provided to accommodate different usage scenarios.

- [**Portable Executable:**](#portable-executable) A ready-to-run, no-install option recommended for quick evaluation and standard usage.
- [**Python Package Installation:**](#python-package-installation) This option installs DDA-BERT via pip and enables users to modify the source code as needed, for example, to support rescoring outputs from additional database search engines. In this mode, DDA-BERT can be used directly as a Python package within custom workflows.
- [**Docker installation**](#docker-installation) Recommended for users who prefer containerized environments or require reproducible deployments.

### Option A. Portable Executable

#### Step 1. Download

Download the latest DDA-BERT portable executable and accompanying test files from the official release page at [https://guomics.com/software/DDA-BERT](https://guomics.com/software/DDA-BERT/downloads.html).

#### Step 2. Run
Unzip the downloaded archive and execute the following command in a terminal:
```shell
cd DDA-BERT; 
./dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output=/out/
```

### Option B. Python Package Installation
It is strongly recommended to install DDA-BERT in an isolated Conda environment.  
⏱️ Estimated setup time: **~10–15 minutes**

**Test data:** demo_data/HeLa_digest_SPME_1ng_1.mzML and demo_data/HeLa_digest_SPME_1ng_1.raw. You can also use your own .raw and .mzML files.

#### Prerequisites

#### Step 1. Download jdk11 from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/DDA-BERT/jdk-11.0.26.zip), unzip and move to project root directory.
> **⚠️Note**: This project was built using **FragPipe v22**, which includes the following core components: **MSFragger v4.1**, **Philosopher v5.1.1**, **diaTracer v1.1.5**, **IonQuant v1.10.27**. If you're using a different version of FragPipe, make sure to download compatible versions of these tools and properly configure your environment to ensure smooth execution of the analysis pipeline.

You can download FragPipe from [here](https://guomics-share.oss-cn-shanghai.aliyuncs.com/SOFTWARE/DDA-BERT/FragPipe22_0.zip) with the following core components.


#### Step 2. Clone the DDA-BERT repository and set up an isolated Conda environment:

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

#### Running DDA-BERT

DDA-BERT provides flexible analysis modes, supporting both an an end-to-end integrated workflow and a modular usage pattern, including:

**Option A. "one-stop" workflow**  
A complete one-stop analysis pipeline encompassing database searching, data preprocessing and cleaning, PSM rescoring, FDR control, and protein inference.
##### Thermo (.raw) data
```shell
dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output=/out/
```

##### Bruker timsTOF (.d) data
```shell
dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output=/out/ --engines=sage --mass-format=d
```

##### Sciex (.wiff) data
Sciex .wiff files should first be converted to .mzML format. Once converted, run the following command:
```shell
dda-bert assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output=/out/ --engines=sage --mass-format=wiff
```

**Option B. Modular rescoring and inference workflow**  
Support for rescoring PSMs from existing database search results, followed by FDR control and protein inference, enabling seamless integration with different search engines or established proteomics workflows.
##### Thermo (.raw) data
```shell
dda-bert score --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --sage-file-dir=xxx --fp-file-dir=xxx --ap-file-dir=xxx --output=/out/
```
##### Bruker timsTOF (.d) data
```bash
dda-bert score --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --sage-file-dir=xxx --output=/out/ --engines=sage --mass-format=d
```
##### Sciex (.wiff) data
```bash
dda-bert score --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --sage-file-dir=xxx --output=/out/ --engines=sage --mass-format=wiff
```
### Option C. Docker Installation

DDA-BERT is available as a self-contained Docker image that includes all required dependencies and runtime environments. This option enables users to run DDA-BERT without manual environment configuration and is well suited for reproducible and portable deployments.

### Prerequisites

* **Docker Engine** installed and running on a Linux system
* Sufficient permissions to pull images from Docker Hub and run containers

Ensure that the Docker service is active before proceeding.

**Step 1. Pull the Docker Image**  
Pull the pre-built DDA-BERT image from Docker Hub:
   ```bash
   docker pull guomics2017/dda-bert:v3.1
   ```
**Step 2. Run the Container**
Docker-based execution uses the same input parameters and argument conventions as described in the `Running DDA-BERT` section above. The specific arguments depend on the mass spectrometry data type (e.g., Thermo .raw, Bruker .d, or converted Sciex .mzML files).
The examples below illustrate how to run the container by mounting a local directory for data access, while using the same input parameters as in the native execution mode.

   ```bash
   docker run --rm -v /home/user/DDA-BERT:/home/test_data guomics2017/dda-bert:v3.1 assess --mzml-paths=/data/example.mzML --fasta=/data/example.fasta --output=/out/
   ```
In this example, the local directory /home/user/DDA-BERT is mounted into the container and used as the working directory for input and output files.
