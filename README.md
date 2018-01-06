# Differential Expression Analysis Workflow        
*By Nikhil Gowda & Sam Chen*

---

### This workflow utilizes the following versions of each tool:
 - FastQC v0.11.5                                           
 - Trimmomatic-0.36                                         
 - STAR-2.5.2b                                              
 - Java (JRE) 1.8.0_131                                     
 - Picard 2.6.0                                             
 - subread-1.5.1                                            
 - R 3.3.1      
Please make sure these are all installed within the scope of the bash environment at runtime.
The program will look for all execpt for java in the bin directory within the project directory.
---

### This workflow utilizes the following reference files:      
 - Reference Genome in fasta format                         
 - Genome Annotations File in gtf format                    
 - Gene Associations File in GAF2.0 format                  
 - RNAseq Data in fastq format  

This workflow expects that the fastq files for different treatment groups are separated in disparate folders.

---

### How to use

#### Run
Run this workflow anytime by running `run.sh` with the bash command. (ie `bash path/to/run.sh`)
The workflow will ask for the location of the required files. The settings can be saved as a `.config` file.

After the first run the workflow can be run with the syntax `bash run.sh group_A_fastq_folder_path group_b_fastq_folder_path`

#### Modifications to the workflow
* Modifications can be made to the parameters provided to workflow tools such as FastQC and STAR by modifying functions.sh
* Modifications to the threshold FDR and FC values can be made in the first few lines of limmavoom.R
