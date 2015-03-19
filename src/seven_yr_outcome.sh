# This tutorial assumes that you have some basic knowledge of working with Linux command-line interface including basic file and folder manipulation operations. It also assumes that you have knowledge of installing packages and compiling source code on Linux as well as setting environment variables for each of the programs used. 

# To proceed ahead with this tutorial, you will need to have the following 5 programs installed/available on a server running Ubuntu 12.04. All analyses will take place via command-line and unfortunately, there is no GUI option available. The URLs provided in front of each program give more description and instructions on downloading/compiling them. 

# QIIME 1.8.0: https://github.com/qiime/qiime-deploy 
# FLASH v1.2.7: http://ccb.jhu.edu/software/FLASH/
# Usearch 64-bit v7.0.1001_i86linux64: http://www.drive5.com/usearch/download.html 
# bcl2fastq v02.14.01.07: http://support.illumina.com/downloads/bcl2fastq-conversion-software-v215.html 
# BIOM : http://biom-format.org/ 


# QIIME consists of a bunch of Python scripts along with various third-party softwares that can all be installed via qiime-deply script (as described on the URL shared above).  
# FLASH is a C program that you would need to compile from source on your Ubuntu machine.
# Usearch is provided as a pre-compiled binary by its developer, so you only need to download it and give it appropriate permissions to be executed. You will need to use a 64-bit copy of Usearch for the analyses.
# bcl2fastq is provided as bundled package by Illumina for machines running CentOS. However, it can be imported into Ubuntu using the alien package. For more details on how to do it, see the link (http://www.howtogeek.com/howto/ubuntu/install-an-rpm-package-on-ubuntu-linux/) 
# BIOM will be installed when QIIME is installed, so no need to install BIOM separately (which is good news). 



# Attach the hard-drive to the server/computer running Ubuntu 12.04. It will show up as URECA_Data in /media folder of Ubuntu 12.04
# Type the following command to go into the folder called "tutorial" where you can run all the analyses steps yourself (through this tutorial). 

cd /media/URECA_Data/7_yr_outcome/tutorial/ 


# You will end up in the main working directory where all the files will be present needed for you analysis
# The raw data produced from Illumina's NextSeq DNA sequencer is present in the directory called "150129_NS500170_0018_AH2HGYAFXX/"
# Go into the directory:

cd 150129_NS500170_0018_AH2HGYAFXX/


# Run bcl2fastq (to convert bcl files to FASTQ files)


bcl2fastq --input-dir Data/Intensities/BaseCalls/ --output-dir Undetermined/ 

# Move the output folder up a directory, go into it and unzip the output files

mv -f Undetermined ../
cd Undetermined/
gunzip *.fastq.gz



# Make directories for data generated via each of the 4 lanes of Illumina's NextSeq in the Undetermined folder

mkdir L001/

mkdir L002/
mkdir L003/
mkdir L004/

# Grant all read/write/execute permissions to the each lane folder

chmod -R 777 L00* 


# Move corresponding lane folders 

mv Undetermined_S0_L001_R* L001/
mv Undetermined_S0_L002_R* L002/
mv Undetermined_S0_L003_R* L003/   
mv Undetermined_S0_L004_R* L004/ 


#Open four separate terminals, with each terminal corresponding to data from each of the 4 lanes of NextSeq:

# LANE 1:

# Change into directory of Lane 1:
cd /media/URECA_Data/7_yr_outcome/tutorial/Undetermined/L001

# Run FLASH program to stitch together Read 1 and Read 2 for 16S rRNA gene:

flash -r 153 -f 253 -s 30 -m 25 -d ./Output_folder_L001/ -q Undetermined_S0_L001_R1_001.fastq Undetermined_S0_L001_R2_001.fastq  

# Move into flash output folder:
cd Output_folder_L001/ 

# Use sed command to get header of FASTQ files as produced by FLASH:

sed -n '1~4'p out.extendedFrags.fastq > FLAShReads.txt

# Go one directory up:

cd ..

# Extract barcodes from FASTQ file using QIIME script called "extract_barcodes.py".  More info on (http://qiime.org/scripts/extract_barcodes.html):

extract_barcodes.py -f Undetermined_S0_L001_R1_001.fastq -c barcode_in_label -s ':' -l 12 -o extracted_barcodes_L001/


# Use a custom python script to order the extracted barcode list the same as those produced by FLASH output

python /media/URECA_Data/7_yr_outcome/tutorial/filter_fasta_keep_order.py extracted_barcodes_L001/barcodes.fastq Output_folder_L001/FLAShReads.txt

# Split libraries to remove barcodes and demultiplex samples using the QIIME script called "split_libraries_fastq.py". More info on (http://qiime.org/scripts/split_libraries_fastq.html):

split_libraries_fastq.py -i Output_folder_L001/out.extendedFrags.fastq -b Index_filtered_ordered.fastq -m /media/URECA_Data/7_yr_outcome/tutorial/7yr_outcome_map.txt -o Split_library_L001/ --rev_comp_barcode --barcode_type 12 -q 30 -s 0



# LANE 2 (Repeat same steps as that of LANE 1. Because input/output file names are different, hence different commands used):
cd /media/URECA_Data/7_yr_outcome/tutorial/Undetermined/L002
flash -r 153 -f 253 -s 30 -m 25 -d ./Output_folder_L002/ -q Undetermined_S0_L002_R1_001.fastq Undetermined_S0_L002_R2_001.fastq 
cd Output_folder_L002/  
sed -n '1~4'p out.extendedFrags.fastq > FLAShReads.txt
cd ..
extract_barcodes.py -f Undetermined_S0_L002_R1_001.fastq -c barcode_in_label -s ':' -l 12 -o extracted_barcodes_L002/
python /media/URECA_Data/7_yr_outcome/tutorial/filter_fasta_keep_order.py extracted_barcodes_L002/barcodes.fastq Output_folder_L002/FLAShReads.txt
split_libraries_fastq.py -i Output_folder_L002/out.extendedFrags.fastq -b Index_filtered_ordered.fastq -m /media/URECA_Data/7_yr_outcome/tutorial/7yr_outcome_map.txt -o Split_library_L002/ --rev_comp_barcode --barcode_type 12 -q 30 -s 30000000



# LANE 3 (Repeat same steps as that of LANE 1. Because input/output file names are different, hence different commands used):
cd /media/URECA_Data/7_yr_outcome/tutorial/Undetermined/L003
flash -r 153 -f 253 -s 30 -m 25 -d ./Output_folder_L003/ -q Undetermined_S0_L003_R1_001.fastq Undetermined_S0_L003_R2_001.fastq  
cd Output_folder_L003/ 
sed -n '1~4'p out.extendedFrags.fastq > FLAShReads.txt
cd ..
extract_barcodes.py -f Undetermined_S0_L003_R1_001.fastq -c barcode_in_label -s ':' -l 12 -o extracted_barcodes_L003/
python /media/URECA_Data/7_yr_outcome/tutorial/filter_fasta_keep_order.py extracted_barcodes_L003/barcodes.fastq Output_folder_L003/FLAShReads.txt
split_libraries_fastq.py -i Output_folder_L003/out.extendedFrags.fastq -b Index_filtered_ordered.fastq -m /media/URECA_Data/7_yr_outcome/tutorial/7yr_outcome_map.txt -o Split_library_L003/ --rev_comp_barcode --barcode_type 12 -q 30 -s 500000000



# LANE 4 (Repeat same steps as that of LANE 1. Because input/output file names are different, hence different commands used):
cd /media/URECA_Data/7_yr_outcome/tutorial/Undetermined/L004
flash -r 153 -f 253 -s 30 -m 25 -d ./Output_folder_L004/ -q Undetermined_S0_L004_R1_001.fastq Undetermined_S0_L004_R2_001.fastq 
cd Output_folder_L004/ 
sed -n '1~4'p out.extendedFrags.fastq > FLAShReads.txt
cd ..
extract_barcodes.py -f Undetermined_S0_L004_R1_001.fastq -c barcode_in_label -s ':' -l 12 -o extracted_barcodes_L004/
python /media/URECA_Data/7_yr_outcome/tutorial/filter_fasta_keep_order.py extracted_barcodes_L004/barcodes.fastq Output_folder_L004/FLAShReads.txt
split_libraries_fastq.py -i Output_folder_L004/out.extendedFrags.fastq -b Index_filtered_ordered.fastq -m /media/URECA_Data/7_yr_outcome/tutorial/7yr_outcome_map.txt -o Split_library_L004/ --rev_comp_barcode --barcode_type 12 -q 30 -s 1000000000



# Once these four terminal steps are over, switch to any one terminal and follow the commands:

cd /media/URECA_Data/7_yr_outcome/tutorial
mkdir all_lanes
chmod 777 all_lanes/ 
cd all_lanes/


mv ../Undetermined/L001/Split_library_L001/seqs.fna seqs_L1.fna 
mv ../Undetermined/L002/Split_library_L002/seqs.fna seqs_L2.fna 
mv ../Undetermined/L003/Split_library_L003/seqs.fna seqs_L3.fna 
mv ../Undetermined/L004/Split_library_L004/seqs.fna seqs_L4.fna 



# Copy "URECA_all_runs_3yr_samples_only_seqs.fna" file to all_lanes directory. This file has sequence data of all individuals (238) for which we had 3 year outcome data.

mv ../URECA_all_runs_3yr_samples_only_seqs.fna .


# Merge all FASTA files into one for downstream analysis:

cat *.fna > URECA_all_runs_7year_seqs.fna


# The FASTA file called "URECA_all_runs_3yr_samples_only_seqs.fna" was created using exactly the same steps used to create the "URECA_all_runs_7year_seqs.fna" but using the L001 to L004 data of samples corresponding to 3 year outcome.

# This section deals with steps used for creation of OTU table from FASTA files using Usearch. We have renamed our Usearch binary as "usearch_64_bit" and run it using ./usearch_64_bit command. It assumes that Usearch is copied in in: /media/URECA_Data/7_yr_outcome/tutorial/all_lanes/Usearch_output. You should also rename your Usearch binary file as "usearch_64_bit", grant it execute permissions and place it in the directory mentioned before, in order to successfully run the commands
# All Usearch commands available at: (http://drive5.com/usearch/manual/uparse_cmds.html )


# Create directory to store output from Usearch


mkdir Usearch_output/
chmod 777 Usearch_output/


# Use a perl script to make the header compatible with Usearch
perl /media/URECA_Data/7_yr_outcome/tutorial/drive5/bmp_uparse_pipeline.pl -i URECA_all_runs_7year_seqs.fna -o URECA_7yr_reads_uparse.fa 


# Use Usearch to dereplicate sequences

./usearch_64_bit -derep_fulllength URECA_7yr_reads_uparse.fa -output URECA_derep_7yrs.fa -sizeout



#Abundance sort and discard singletons
#This sorts your dereplicated reads based on the number of times each was found.  If you want to keep singletons there is a flag for that.
./usearch_64_bit -sortbysize URECA_derep_7yrs.fa -output URECA_derep_7yrs_sorted.fa -minsize 2 


#OTU clustering using UPARSE method
./usearch_64_bit -cluster_otus URECA_derep_7yrs_sorted.fa -otus URECA_otus_picked_7yrs.fa


#Label OTU sequences OTU_1, OTU_2...
python /media/URECA_Data/7_yr_outcome/tutorial/drive5/fasta_number.py URECA_otus_picked_7yrs.fa OTU_>URECA_otus_nochimera_labeled_7yrs.fa 


#Map reads (including singletons) back to OTUs
./usearch_64_bit -usearch_global URECA_7yr_reads_uparse.fa -db URECA_otus_nochimera_labeled_7yrs.fa -strand plus -id 0.97 -uc map.uc


#Create OTU table
python /media/URECA_Data/7_yr_outcome/tutorial/drive5/uc2otutab.py map.uc > URECA_otu_7yr_final_table.txt




# Go back to tutorial directory, create a new directory called "otu_table_finalization", change permissions to it and go into the new directory:


cd /media/URECA_Data/7_yr_outcome/tutorial/
mkdir otu_table_finalization/

chmod 777 otu_table_finalization/
cd otu_table_finalization/


# Convert OTU table generated via Usearch into BIOM format. For more info (http://biom-format.org/ ): 

biom convert -i /media/URECA_Data/7_yr_outcome/tutorial/all_lanes/Usearch_output/URECA_otu_7yr_final_table.txt -o otu_table_uparse.biom --table-type="otu table"


# Assign taxonomy to OTUs using uclust method using QIIME script "assign_taxonomy.py" (use the file “URECA_otus_nochimera_labeled_7yrs.fa” from Usearch as input file). For more info (http://qiime.org/scripts/assign_taxonomy.html):

assign_taxonomy.py -i /media/URECA_Data/7_yr_outcome/tutorial/all_lanes/Usearch_output/URECA_otus_nochimera_labeled_7yrs.fa -o output_otu_taxonomy -r /media/URECA_Data/7_yr_outcome/tutorial/Greengenes_Database_May_2013/97_otus.fasta -t /media/URECA_Data/7_yr_outcome/tutorial/Greengenes_Database_May_2013/taxonomy/97_otu_taxonomy.txt



# Add metadata (taxonomy) to OTU table. For more info (http://biom-format.org/ ): 


biom add-metadata -i otu_table_uparse.biom -o otu_table_with_tax_7yr_outcome.biom --observation-metadata-fp output_otu_taxonomy/URECA_otus_nochimera_labeled_7yrs_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence


# Summarize OTU table.  For more info (http://biom-format.org/ ): 

biom summarize-table -i otu_table_with_tax_7yr_outcome.biom -o biom.summary


# Filter the non-corresponding samples from OTU table to keep only the relevant 279 samples (238+41). For more info (http://qiime.org/scripts/filter_samples_from_otu_table.html):

filter_samples_from_otu_table.py -i otu_table_with_tax_7yr_outcome.biom -o filtered_279_samples_with_tax_7yr_outcome.biom --sample_id_fp /media/URECA_Data/7_yr_outcome/Analysis_After_Usearch_021315/all_279_sample_ids_for_filtering.txt


# Summarize filtered OTU table. For more info (http://biom-format.org/ ): 

biom summarize-table -i filtered_279_samples_with_tax_7yr_outcome.biom -o biom_279_samples.summary 


# Discard the two non-amplified samples from OTU table to keep leaving 277 samples in total.  For more info (http://qiime.org/scripts/filter_samples_from_otu_table.html)

filter_samples_from_otu_table.py -i filtered_279_samples_with_tax_7yr_outcome.biom -o final_277_samples_with_tax_7yr_outcome.biom --sample_id_fp /media/URECA_Data/7_yr_outcome/Analysis_After_Usearch_021315/all_277_sample_ids_for_filtering.txt



# BIOM convert for multiply rarefying step.  For more info (http://biom-format.org/ ): 

biom convert -b -i final_277_samples_with_tax_7yr_outcome.biom -o final_277_samples_with_tax_7yr_outcome.txt --header-key taxonomy


# Multiple rarefying at depth of 173841. This is an Rscript written to perform multiple rarefying. More details on this script available at (https://github.com/alifar76/MicroNorm)

Rscript multiple_rarefying_script.R final_277_samples_with_tax_7yr_outcome.txt multiply_rarefied_277_samples_173841.txt 100 canberra 


# The table produced "multiply_rarefied_277_samples_173841.txt" is the locked down, final OTU table table. This OTU table is present in the directory:  /media/URECA_Data/7_yr_outcome/Analysis_After_Usearch_021315
