# Mammalian-evolutionary-rates
This is the project involving the analyses of mammalian genomes to understand the dynamics, causes and impacts of mammalian evolutionary rates 
The intermediate files are presented in the excel table.

Two python scripts are also presented which do not require any major dependency.

1. Read cleaner cleans the read, and can be simply run with the following command.

read_cleaner.py -i Capybara_raw.1.fastq -j  Capybara_raw.2.fastq -1  Capybara_filtered.1.fastq -2  / 
Capybara_filtered.1.fastq -u Capybara_filtered_unpaired.1.fastq –v / Capybara_filtered_unpaired.2.fastq -q 20 -l 50

read_cleaner.py (without any option) or read_cleaner.py -h will give the help for the usage


1. Depth masker masks out unreliable genomic positions.

depth_masker.py -v Capybara_20Q_50bp_d3D30.vcf -o Capybara_20Q_50bp_d3D30.fa  -d 3 /
 –D 30 -Q 30 

