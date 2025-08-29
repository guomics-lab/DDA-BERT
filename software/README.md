
<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://github.com/zhiyuajun/DDA-BERT/blob/main/DDA-BERT.png" width="240" height="240">
</p>

## Description
DDA-BERT: a computational platform designed for data-dependent acquisition (DDA) proteomics analysis.

## Software
The software and manual can be downloaded from the website https://guomics.com/software/DDA-BERT.
On Linux, download the file from the release. DDA-BERT runs install-free and requires no additional configuration of the environment. 

## Installation
If you want to use DDA-BERT by source code, you can install python and install requirements package.

### Prerequisites
Please make sure you have a valid installation of conda or miniconda. We recommend setting up miniconda as described on their website.

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
•	Operating System: Supports both Windows and Linux operating systems.

•	Processor: A dual-core processor is recommended, but it can run on a single-core processor.

•	Memory: 40GB or more is recommended. If the mass spectrometry files or library files to be identified are large, it is advised to use more memory.

•	Storage: At least 100GB of available hard disk space is recommended.

•	Graphics Card: A 40GB NVIDIA GPU with CUDA support or a V100 32GB GPU is recommended.

## License
This software is licensed under a custom license that allows academic use but prohibits commercial use. For more details, see the LICENSE file.

## Contact
For any questions or licensing inquiries, please contact:
Dr Guo
E-mail: guotiannan@westlake.edu.cn
www.guomics.com
