## Description
DDA-BERT: a computational platform designed for data-dependent acquisition (DDA) proteomics analysis.

## Software
The software and manual can be downloaded from the website https://guomics.com/software/DDA-BERT.
On Linux, download the file from the release. DDA-BERT runs install-free and requires no additional configuration of the environment. 

## Installation
To use DDA-BERT from source, install Python and the required dependencies.

### Prerequisites
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
pip install -r requirements.txt
```

**You need install torch from pytorch (https://pytorch.org/).** It is advisable to install the entire pytorch package and follow the official installation method provided by pytorch.

Specifically, first select the CUDA version according to your own operating system, and then, based on the CUDA version, choose the corresponding installation command to execute. For example, run "pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124".

Linux command-line run
```shell
python main_linux.py --mzml_paths=/data/example.mzML --fasta=/data/example.fasta --output_path=/out/
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
