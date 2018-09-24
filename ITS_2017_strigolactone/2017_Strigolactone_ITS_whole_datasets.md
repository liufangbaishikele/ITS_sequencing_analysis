## ITS analysis using 2017 strigolactone whole dataset

1. Transfer data from local to staton server using ``scp``

```
scp -r -P 22 Local_path_on_PC fliu21@Staton.XX.XX:pathway_for_raw_data_File_inside_staton_server
```

2. Unzip files and re-arrange raw fastq file to be ready for analysis


```
cd /staton/projects/soybean_rhizosphere/ITS_2017_strigolactone/raw_data
for file in *; do echo $file; cd $file; gunzip *; cd ../; done
mkdir ../raw_fastq
for file in *; do echo $file; cd $file; mv * ../../raw_fastq; cd ../; done
```

3. Generate ITS.file and ITS.oligo files in raw sequence file

```
cd ../raw_fastq
ls *R1_001.fastq > R1_list
ls *R2_001.fastq > R2_list
sed 's/\(^.*\)_S.*/\1/g' R1_list > sample-ID
sed 's/-/_/g' sample-ID > sample_ID
vim sample_ID and change G1_ to D14_, G10_ to Max1_ and G14_ to Max2_.
paste -d "\t" sample_ID R1_list R2_list > ITS.file
cp ../mothur_two_samples/rawdata/ITS.oligo ./
```

4. Move UNITE reference to analysis file


5. Double check ITS.batch file for batch based mothur job
