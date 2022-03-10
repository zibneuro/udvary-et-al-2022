### Analysis routines
## Contents
- visualization: scripts to reproduce the figures and tables (Matlab)
- output: output directory for figures and tables
- preprocessing: data preprocessing and format conversion scripts (Matlab, C++)

## System requirements
- Matlab R2020b (recommended)

## Reproduce figures and tables
### Download required data
- Download 'analyzed data' from [BarreCortexInSilico](https://cortexinsilico.zib.de/download) 
- Extract the zip file and copy the folder 'data' to the location of the placeholder file analyzed_data.txt in this folder.

### Run Matlab scripts
- Add folder 'visualization' to Matlab path. Use: 'Set Path', 'Add with Subfolders'.
- Run scripts in folder 'visualization' independently. 
- Inside each script, adjust the 'matlabPath' (default: 'D:\udvary-et-al-2022\analysis\') according to the absolute path of this repository.

## Running preprocessing scripts
### Download required data
- Download 'intermediate data' from [BarreCortexInSilico](https://cortexinsilico.zib.de/download) 
- Extract the zip file and copy the folder 'data' to the location of the placeholder file intermediate_data.txt in folder 'preprocessing'.

### Run Matlab scripts
- Add folder 'preprocessing' to Matlab path. Use: 'Set Path', 'Add with Subfolders'.
- Run scripts in folder 'preprocessing' independently. 
