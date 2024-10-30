# DDA-BERT
# Description
DDA-BERT is a computational software for data-dependent acquisition (DDA) proteomics analysis. It offers a simple, installation-free graphical user interface (GUI) and can also run in command-line mode. Compatible with both Windows and Linux, DDA-BERT utilizes the ultra-fast Sage engine, a specialized proteomics search engine for peptide identification.   
DDA-BERT rescores peptide identifications using a pre-trained Transformer-based model trained on 3,706 DDA files and over 95 million peptide-spectrum matches (PSMs). Currently, the platform supports .mzML and Bruker .d formats. Additional formats, such as Sciex .wiff and Thermo .raw, can be converted to .mzML via the MSConvertGUI tool from ProteoWizard. Future updates will extend direct support to these and other formats.

# Installation
On Windows systems, download and unzip the zip file. Click on dda-bert.exe to run without the installation of any runtimes (e.g., Java, Python) or external dependencies.  
On Linux, download the file from the release. DDA-BERT runs install-free and requires no additional configuration of the environment. 

Hardware Requirements:  
•	Operating System: Supports both Windows and Linux operating systems.  
•	Processor: A dual-core processor is recommended, although it can run on a single-core processor.  
•	Memory: 40GB or more is recommended. For larger mass spectrometry or FASTA files, additional memory is advisable.  
•	Storage: A minimum of 100GB of available hard disk space is recommended.  
•	Graphics Card: A 40GB NVIDIA GPU with CUDA support is recommended.

DDA-BERT executables are available via https://dda-bert.guomics.com/.

## Run the analysis from the command line
1.	Windows  
cd gpt-desktop-win32-x64;
./dda-bert-start.exe --mzml_paths=/data/example.mzML --fasta=/data/example.fasta --output_path=/out/

2.	Linux  
cd dda-bert-start;
chmod -R a+x sage/linux*;
./dda-bert-start --mzml_paths=/data/ example.mzML --fasta=/data/example.fasta --output_path=/out/

## Results
Results are output in CSV format as a comprehensive summary table that is easy to manipulate and interpret, facilitating further biological insights and downstream applications.

# License
This software is licensed under a custom license that allows personal use but prohibits commercial use. For more details, see the LICENSE file.

# Contact
For any questions or licensing inquiries, please contact: Dr Guo E-mail: guotiannan@westlake.edu.cn
www.guomics.com
