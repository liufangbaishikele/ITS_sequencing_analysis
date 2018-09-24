##                               ITS sequencing data analysis using mothur

---
Data set comes from 2017 strigolactone project, with RNAi, Over-expression and GR24 restore treatments
---

**Target region**: ITS2 by ITS3_F and ITS4_R, for arrangement and primer choice of ITS1 and ITS2 please refer [PacBio metabarcoding of Fungi and other eukaryotes: errors,biases and perspectives](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.14776) for more information.

## Analysis process

1) Transfer data to staton server
2) Download UNITE reference follower [UNITE souce](https://unite.ut.ee/repository.php)
3) Come up with reasonable strategy for this data analysis.

    Refer to these links for inspiring comments about ITS analysis pipeline on mothur
    
    i) [Options to run fungi dataset into Mothur #334](https://github.com/mothur/mothur/issues/334)
    
    ii) [How to conduct MOTHUR analysis of fungal ITS sequence data ?](https://www.researchgate.net/post/How_to_conduct_MOTHUR_analysis_of_fungal_ITS_sequence_data#view=5b8b808584a7c142a218a9b0)
  
    iii) [Mirza Abid post about UNITE and ITS pipeline](https://www.researchgate.net/post/Mothur-formatted_UNITE_database)
    
4) Pre-process of sequences before clustering

```
make.contigs(inputdir=/staton/projects/soybean_rhizosphere/ITS_2017_strigolactone/mothur_two_samples/rawdata,outputdir=/staton/projects/soybean_rhizosphere/ITS_2017_strigolactone/mothur_two_samples/analysis,file=ITS.file,oligos=ITS.oligo,processors=12)

summary.seqs(inputdir=/staton/projects/soybean_rhizosphere/ITS_2017_strigolactone/mothur_two_samples/analysis,fasta=ITS.trim.contigs.fasta,processors=12)

screen.seqs(fasta=ITS.trim.contigs.fasta,group=ITS.contigs.groups,maxambig=0,maxlength=428,processors=12)

summary.seqs(fasta=ITS.trim.contigs.good.fasta,processors=12)

unique.seqs(fasta=ITS.trim.contigs.good.fasta)

count.seqs(name=ITS.trim.contigs.good.names,group=ITS.contigs.good.groups)

summary.seqs(count=ITS.trim.contigs.good.count_table,processors=12)

```

**Skip alignment related process**, althougth alignment based clustering will be more accurate but for ITS amplicon, this is not feasible as too much variables of ITS regions between taxons

```
#align.seqs(fasta=ITS.trim.contigs.good.unique.fasta,reference=UNITEv6_sh_97.fasta,processors=12)

#summary.seqs(fasta=ITS.trim.contigs.good.unique.align,count=ITS.trim.contigs.good.count_table,processors=12)

#screen.seqs(fasta=ITS.trim.contigs.good.unique.align,count=ITS.trim.contigs.good.count_table,summary=ITS.trim.contigs.good.unique.summary,start=2,end=17012,maxhomop=8,processors=12)

#summary.seqs(fasta=ITS.trim.contigs.good.unique.good.align,count=ITS.trim.contigs.good.good.count_table,processors=12)

#filter.seqs(fasta=ITS.trim.contigs.good.unique.good.align,vertical=T,trump=.,processors=12)

#unique.seqs(fasta=ITS.trim.contigs.good.unique.good.filter.fasta,count=ITS.trim.contigs.good.good.count_table,processors=8)
processors is not a valid parameter.

#unique.seqs(fasta=ITS.trim.contigs.good.unique.good.filter.fasta,count=ITS.trim.contigs.good.good.count_table)
```
**Continue with pre.cluster**

```
pre.cluster(fasta=ITS.trim.contigs.good.unique.fasta,count=ITS.trim.contigs.good.count_table,diffs=3,processors=12)

chimera.vsearch(fasta=ITS.trim.contigs.good.unique.precluster.fasta,count=ITS.trim.contigs.good.unique.precluster.count_table,dereplicate=t,processors=12)

remove.seqs(fasta=ITS.trim.contigs.good.unique.precluster.fasta,accnos=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta=ITS.trim.contigs.good.unique.precluster.pick.fasta,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,processors=12)


classify.seqs(inputdir=/staton/projects/soybean_rhizosphere/ITS_2017_strigolactone/mothur_two_samples/analysis/pairwise.seq_and_clustering,fasta=ITS.trim.contigs.good.unique.precluster.pick.fasta,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,reference=UNITEv6_sh_97.fasta,taxonomy=UNITEv6_sh_97.tax,cutoff=80,processors=12)

remove.lineage(fasta=ITS.trim.contigs.good.unique.precluster.pick.fasta,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table,taxonomy=ITS.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.taxonomy,taxon=unknown)

summary.tax(taxonomy=ITS.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table)

```

##  Here we need to bypass the alignment based clustering and alternatively 


#### 4-1) ``pairwise.seqs`` to cluster.seqs as the former command will generate a column fomatted distance matrix. Then, the ``cluster`` will do the clustering to generate OTU table.

```
pairwise.seqs(fasta=ITS.trim.contigs.good.unique.precluster.pick.pick.fasta,cutoff=0.1,align=needleman,output=lt,countends=T,calc=onegap,processors=12)

#--- Clustering using opticlust algrithm

cluster(phylip=ITS.trim.contigs.good.unique.precluster.pick.pick.phylip.dist,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,method=opti,cutoff=0.03,processors=24)

#---Using generated list file to create OTU table

make.shared(list=ITS.trim.contigs.good.unique.precluster.pick.pick.phylip.opti_mcc.list,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table)

classify.otu(list=ITS.trim.contigs.good.unique.precluster.pick.pick.phylip.opti_mcc.list,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=ITS.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,label=0.03)

```

#### 4-2) Vsearch greedy clustering based OTU clustering, using ``cluster`` with ``dgc`` (distance based greedy clustering) method

Run cluster.seq command with vsearch based clustering mthod, e.g., agc (abundance based greedy clusterin) and dgc (distance based greedy clustering, with the dgc as the default)

a) First, copy fasta, count table and taxonomy file to vsearch_clustering directory. Here the fasta is after chimera screen and unknown fungi removing. classify.seqs were done and all unknown sequence were removed from the fasta file

b) Doing clustering first, here the cutoff need to set up to 0.03 as the generated list file will be OTU with 97% similarity. Here the method is using 

```
cluster(fasta=ITS.trim.contigs.good.unique.precluster.pick.pick.fasta,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,method=dgc,cutoff=0.03,processors=48)

make.shared(list=ITS.trim.contigs.good.unique.precluster.pick.pick.dgc.list,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,label=0.03)

classify.otu(list=ITS.trim.contigs.good.unique.precluster.pick.pick.dgc.list,count=ITS.trim.contigs.good.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=ITS.trim.contigs.good.unique.precluster.pick.UNITEv6_sh_97.wang.pick.taxonomy,label=0.03)
```

