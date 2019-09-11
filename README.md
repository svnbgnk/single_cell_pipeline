## Single-Cell-Nanopore – nanopore read assignment

This modified version of the program Flexbar assigns barodes to nano nanopore sequencing data by using the
    efficient exact overlap alignments using SIMD and multicore parallelism provided by flexbar and SeQan. It removes the
    provided adapter sequence and assign barcodes provided in a separate fasta, which could be determine by the
    chromium pipeline.

[comment]: # Mention seqan

[comment]: # Refer to the [manual](https://github.com/seqan/flexbar/wiki) or contact [](https://github.com/jtroehr) for support with this application.


### References
[comment]: #  Add additional references
Johannes T. Roehr, Christoph Dieterich, Knut Reinert:  
Flexbar 3.0 – SIMD and multicore parallelization. Bioinformatics 2017.

See article on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28541403)

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich:  
Flexbar – flexible barcode and adapter processing for next-generation sequencing platforms. Biology 2012.

See article on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24832523)

Cloning:
	git clone https://github.com/svnbgnk/single_cell_pipeline.git

Use these commands for building:
	mkdir flexbar-build
	cd flexbar-build
	cmake ../flexbar
	make


### Binaries

For execution of provided Flexbar binaries, the corresponding TBB library has to be available. Downloads contain the library file for runtime. Follow the platform specific instructions below.

### Program usage

SingleCellPipe needs a bam alignment file, a barcode whitelist, gtf file containing regions of interest and a primer sequence. Additionally, the target name and further options can be specified. 

[comment]: #	singleCellPipe -r reads -w [-b barcodes] [-a adapters] [options]
	singleCellPipe -r nanoporeAlignmentFile -w CellbarcodeWhitelist -rf regionsOfInterest -t Target -as "Primer Sequence" [options]


Refer to the help screen `singleCellPipe -h` for more information.

### Advanced options

singleCellPipe -r nanoporeAlignmentFile -w CellbarcodeWhitelist -rf regionsOfInterest -t 
 -as "Primer Sequence" -n 40 -N 3 -fq -rm


### Input Examples

##Cellbarcode Whitelist
Ids will be used in log file to identify the align cellbarcode.
```
>CellbarcodeID1
CGTCCATCAGTGCGCT
>CellbarcodeID2
AATCGTGCATCATTTC
>CellbarcodeID3
AGTGACTCAGTGCGCT

```


##Regions file
For the regions file only the first, fourth and fifth column are considered. 
```
1	ensembl_havana	gene	65419	71585	.	+	.	gene_id "ENSG00000186092"; gene_version "5"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
1	ensembl_havana	gene	450703	451697	.	-	.	gene_id "ENSG00000284733"; gene_version "1"; gene_name "OR4F29"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
```
##Primer Sequence
The primer sequence is given as a plain string in qutations:
```
 -as "CTACACGACGCTCTTCCGATCT"
```

###Output 
The Cellbarcode CellbarcodeID1 has 3 valid alignments to barcodes of the whitelist. The difference two the second best score is 7. 
This is summarized with Best alignment (3/7). 

```
Alternative alignment:
Sequence removal: left side
  query id         CellbarcodeID2
  query pos        0-16
  read id          63808c32-5d28-4996-bc4d-7d877753a536_end1_Flexbar_removal_cmdline_Flexbar_removal_AATCGTGCATCATTTC
  read pos         3-19
  score            7
...


Alternative alignment:
Sequence removal: left side
  query id         CellbarcodeID3
  query pos        0-17
  read id          63808c32-5d28-4996-bc4d-7d877753a536_end1_Flexbar_removal_cmdline_Flexbar_removal_AGTGACTCAGTGCGCT
  read pos         0-17
  score            9
...


Best alignment (3/7):
Sequence removal: left side
  query id         CellbarcodeID1
  query pos        0-16
  read id          63808c32-5d28-4996-bc4d-7d877753a536_end1_Flexbar_removal_cmdline_Flexbar_removal_CGTCCATCAGTGCGCT
  read pos         0-16
  score            16
  overlap          16
  errors           0
  error threshold  3.84
  remaining read   TTTGGTTAG
  remaining qual   =<-3;/=/*
  barcode seq      CGTCCATCAGTGCGCT
  barcode qual     ('-.>99;2(8*+,-7


  Alignment:

      0     .    :    .  
        CGTCCATCAGTGCGCT
        ||||||||||||||||
        CGTCCATCAGTGCGCT
```
