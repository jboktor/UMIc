# UMIc

UMIc is a framework, written in [R](https://www.r-project.org/), implementing a new proposed method for UMI deduplication and reads correction. The method works at nucleotide level, taking into account the frequency and the mean quality of each base. 

## Getting started

### Prerequisites

The packages needed to be installed, in order to run the project are:
- from CRAN

```
install.packages(c("tidyverse", "data.table", "stringr", "stringdist"))
```

- from Bioconductor
```
BiocManager::install(c("Biostrings", "ShortRead"))
```

### Installing

The project can be downloaded using git:

```
git clone https://github.com/BiodataAnalysisGroup/UMIc
```

### Running the project

The framework consists of three scripts:
- ```UMIsProject.R``` 
- ```casesWorkflows.R```
- ```functions.R```

In order to run the project, set the [R](https://github.com/BiodataAnalysisGroup/UMIc/tree/master/R) folder as your working directory, set the parameters input parameters in the main script ```UMIsProject.R``` and then use the following command:
```
source("UMIsProject.R")
```

The project provides example input datasets and their outputs, for testing purposes. The folder [data](https://github.com/BiodataAnalysisGroup/UMIc/tree/master/data) includes example datasets for all three scenarios in their corresponding subfolders. Each subfolder icludes the fastq files and a Readme.md file with the parameter values, used to generate the files in folder [outputs](https://github.com/BiodataAnalysisGroup/UMIc/tree/master/outputs). The user must provide the input and output folders' filepaths. 

### Inputs
Before running the project, the user must set the appropriate input parameters in the main script ```UMIsProject.R```.

The later has the following inputs:
- ```pairedData```: boolean variable that indicates, whether data are paired ```T``` or single ```F```. 
- ```UMIlocation```: variable that indicates, whether UMI is located only in Read1 ```R1``` or Read1 and Read2 ```R1 & R2```.
- ```UMIlength```: the length of the UMI sequence.
- ```sequenceLength```: the length of the read sequence. 
- ```countsCutoff```: min read counts per UMI, for initial data cleaning.
- ```UMIdistance```: max UMI distance for UMI merging.
- ```sequenceDistance```: max sequence distance for UMI merging.
- ```inputsFolder```: name or filepath of the inputs folder.
- ```outputsFolder```: name or filepath of the outputs folder, default value is ```UMIc_output```.
 
The input data must be provided in fastq files and it is assumed that the UMI is placed at the beginning of each sequence. The library preparation step of the input files must be genarated using the same protocol and fulfil the same input parameters described above.


### Outputs 
The output data are stored also in fastq files, named the same as the input files with an added ```_corrected``` suffix and the name of the folder can be provided by the user. The files contain the corrected sequences (without the UMI) and their quality. It is worth mentioning that the new sequence ID is constructed by combining the ID of one of the input sequences, that has that same UMI, and the UMI itself.

The framework also produces a csv file with all the information of the output fastq files and extra information, that can help return from the output sequences to their corresponding input sequences. The file is named the same as the Read1 fastq file with an added ```_summary_table``` suffix.

For more details, please refer to the [wiki](https://github.com/BiodataAnalysisGroup/UMIc/wiki).


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
