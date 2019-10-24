# On-boarding

Hello and welcome!  
This document is a collection of the tools and tasks that the lab is interested in.

# Gene expression
## Differential expression, co-expression and differential co-expression
<img src="./imgs/schematic.png" width="1234" height="314" title="schematic"> 

![summary](imgs/schematic.png "schematic")

![summary](imgs/netagg.png "aggre")

![perf1](imgs/assess.png "egad") ![perf2](imgs/coexpp.png "performance")



# Databases and repositories  
## Gene expression data 
Databases for gene expression data in raw and parsed form. 

[GEO](https://www.ncbi.nlm.nih.gov/geo/) 

[SRA](https://www.ncbi.nlm.nih.gov/sra)

[ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)

[ENA](https://www.ebi.ac.uk/ena)

[Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home)

[Expression Atlas](https://www.ebi.ac.uk/gxa/home)

[Allen Brain Atlas](https://portal.brain-map.org/)

## Annotations
[GENCODE](https://www.gencodegenes.org/)

[Human](https://www.gencodegenes.org/human/)

[Mouse](https://www.gencodegenes.org/mouse/)

[modENCODE (for fly and worm)](http://www.modencode.org/) 

[ENSEMBL](http://useast.ensembl.org/index.html)  

[Other/all species](http://useast.ensembl.org/info/about/species.html)

 
## Core datasets 
[GTEx](https://gtexportal.org/home/)

[GEUVADIS](https://www.ebi.ac.uk/Tools/geuvadis-das/)

[ENCODE](https://www.encodeproject.org/) 

[BrainSpan](https://www.brainspan.org/) 

## Processed expression data 
[Recount2](https://jhubiostatistics.shinyapps.io/recount/) 

[GEMMA](https://gemma.msl.ubc.ca/)

[ARCHS4](https://amp.pharm.mssm.edu/archs4/)

## Co-expression databases 
[COEXPRESdb](https://coxpresdb.jp/)

[HumanBase](https://hb.flatironinstitute.org/)


## Ontologies
Ontologies have vocabularies (term -> term) and annotation (gene -> term) relationships. Some useful and key ontologies:  
[Gene Ontology](http://geneontology.org/) 
Related [work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389513/). 

[Human phenotype ontology (HPO)](https://hpo.jax.org/app/)

[Cell ontologies](http://www.obofoundry.org/ontology/cl.html)

[Experimental factor (EFO)](https://www.ebi.ac.uk/efo/)

Others (and data) can be found here at the [Harmonizome](https://amp.pharm.mssm.edu/Harmonizome/). 
The [OBO foundry](https://github.com/OBOFoundry/purl.obolibrary.org/) has the standard vocabularies. 

## Pathways 
[KEGG](https://www.genome.jp/kegg/)

[Reactome](https://reactome.org/)

[Biocarta](https://amp.pharm.mssm.edu/Harmonizome/dataset/Biocarta+Pathways). The original site seems to be dead/down. 


## Protein data 
[biogrid](https://thebiogrid.org/) 

[STRING](https://string-db.org/)

[I2D](http://ophid.utoronto.ca/ophidv2.204/) 

[HIPPIE](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/)

[Interpro](https://www.ebi.ac.uk/interpro/)

[HuRI (human reference or CCSB)](http://interactome.baderlab.org/download)

[Human Protein Atlas](https://www.proteinatlas.org/)

## Genomic data
[1000 genomes](https://www.internationalgenome.org/)

[ExAC](http://exac.broadinstitute.org/)

[gnoMAD](https://gnomad.broadinstitute.org/)

[TOPmed](https://www.nhlbiwgs.org/) which can also be accessed [here](https://bravo.sph.umich.edu/freeze5/hg38/). 

[dbGAP](https://www.ncbi.nlm.nih.gov/gap/)

### Key tools
[GATK](https://software.broadinstitute.org/gatk/)

[Samtools](http://www.htslib.org/) 

[Tabix](http://www.htslib.org/doc/tabix.html)

[IGVtools](https://software.broadinstitute.org/software/igv/igvtools)


# Gene expression analysis
## Microarray
Notes [here](https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor) 

## RNA-sequencing 
Some useful notes [here](https://bioinformatics-core-shared-training.github.io/RNAseq-R/) and [here](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)

### Bulk 
List of tools [here](https://en.wikipedia.org/wiki/List_of_RNA-Seq_bioinformatics_tools) 
And some others like [fastX](http://hannonlab.cshl.edu/fastx_toolkit/) and [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are good for QC. 

### Single-cell 
List of all tools [here](https://www.scrna-tools.org). Some key tools include [Seurat](https://satijalab.org/seurat/) and the [Hemberg lab's course](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html).

### Other 
[Biojupies](https://amp.pharm.mssm.edu/biojupies/). 

# Gene set enrichment 
## Key tools 
[ermineJ](https://erminej.msl.ubc.ca/) 

[GSEA](http://software.broadinstitute.org/gsea/index.jsp)

[DAVID](https://david.ncifcrf.gov/) 

[GEO2Enrichr](https://amp.pharm.mssm.edu/g2e/)

# Alignment tools 
## STAR
Github [here](https://github.com/alexdobin/STAR) and [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). 
Reference [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) and [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/).

## Kallisto
Github [here](https://pachterlab.github.io/kallisto/) and [tutorial](https://pachterlab.github.io/kallisto/starting). 
Reference [here](https://www.nature.com/articles/nbt.3519). 

## Salmon
Github [here](https://combine-lab.github.io/salmon/) and [manual](https://combine-lab.github.io/salmon/getting_started/). 
Reference [here](https://www.nature.com/articles/nmeth.4197). The single cell version (Alevin) can be found [here](https://salmon.readthedocs.io/en/latest/alevin.html) and [ref](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y). 

## Bowtie2
[Source](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). 
References [here](https://academic.oup.com/bioinformatics/article/35/3/421/5055585), [here](https://www.nature.com/articles/nmeth.1923) and [here](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25). 


# Model organisms 
## Orthology 
[Homologene](https://www.ncbi.nlm.nih.gov/homologene)

[OrthoDB](https://www.orthodb.org/)

[BUSCO genes](https://busco.ezlab.org/) 

## Species of interest  
- Mouse (Mus musculus, 10090) at [JAX](http://www.informatics.jax.org/)
- Yeast	(Saccharomyces cerevisiae,4932) at [yeastgenome](https://www.yeastgenome.org/) or (Schizosaccharomyces pombe, 284812) at [pombase](https://www.pombase.org/)
- Fly	(Drosophila melanogaster,7227)	at [flybase](https://flybase.org/) 
- Maize	(Zea mays, 4577) at [maizegdb](https://www.maizegdb.org/), [gramene](http://www.gramene.org/), or at [ensembl](https://plants.ensembl.org/Zea_mays/Info/Index) 
- Arabidopsis	(Arabidopsis thaliana, 3702) [here](https://www.arabidopsis.org/) or [AtGDB](http://www.plantgdb.org/AtGDB/)
- Worm (Caenorhabditis elegans,6239) at [wormbase](https://www.wormbase.org/)
- Zebrafish (Danio rerio, 7955) at [zfin](https://zfin.org/)
- Frog (Xenopus laevis, 8355) at [xenbase](http://www.xenbase.org/entry/)
- Armadillo	(Dasypus novemcinctus,9361) at [here](http://www.xenarthrans.org/) and [ensembl]().
- Naked mole rats	(Heterocephalus glaber, 10181) at [here](http://www.naked-mole-rat.org/)


# Lab bits and pieces  
https://gillislab.github.io/publications/

https://gillislab.github.io/background/ 



