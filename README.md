# p-dindel

Author: Mohammad Shabbir Hasan
PhD Student, Department of Computer Science
Virginia Tech, Blacksburg, VA 24060, USA.
E-mail: shabbir5@vt.edu
===========================================================

## System Requirement

-GCC and Python.
-P-Dindel has been tested with GCC version 4.6.3 and Python version 2.7.3 in the Ubuntu 12.04 LTS operating system.

## Compiling

-To compile p-dindel.cpp use the following command
  g++ p-dindel.cpp -std=c++0x -lpthread -lrt -o p-dindel

-To compile p-dindel_pooled.cpp use the following command
  g++ p-dindel_pooled.cpp -std=c++0x -lpthread -lrt -o p-dindel_pooled
  
## Execution of the program

-To run P-Dindel for diploid sample, you need to use the following command:
  ./p-dindel -i input_bam_file_name -o output_file_name -r reference_fasta_file_name

-To run P-Dindel for pooled samples, you need to use the following command:
  ./p-dindel_pooled -i input.txt -o output_file_name -r reference_fasta_file_name
  
  Here input.txt contains the list of bam file for each sample in the pool (1 bam file in each line).
