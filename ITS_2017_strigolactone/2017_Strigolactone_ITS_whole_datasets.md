## ITS analysis using 2017 strigolactone whole dataset

1. Transfer data from local to staton server using ``scp``

```
scp -r -P 22 Local_path_on_PC fliu21@Staton.XX.XX:pathway_for_raw_data_File_inside_staton_server
```

2. Unzip files and re-arrange raw fastq file to be ready for analysis

```
for file in *; do echo $file; cd $file; gunzip *; cd ../; done
```

3. Generate ITS.file and ITS.oligo files in raw sequence file

```

```

4. Move UNITE reference to analysis file
5. Double check ITS.batch file for batch based mothur job
