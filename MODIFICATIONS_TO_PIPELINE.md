# Validation pipeline
## Overview and background
The idea of the current pipeline is to “validate” eDNA signals that we get by classifying (Kraken2) a huge metagenomics dataset spanning 516 samples (~150-200 million PE reads per sample) across 10 years. Our main hypothesis is that a true positive species should generate reads randomly and uniformly from across its genome. For any given species in our dataset, we want to use assumption by mapping the reads classified to it to its reference genome and then calculating the correlation between the length of a contig and the number of reads mapping to it. The correlation should be very high for a true positive species. However, there are complicating factors. For example, the ability of kraken to classify the reads to a species depends on how many minimizers are associated with that taxon, and the distribution of those. If a species has many sibling species (within the same genus), it will share much of the reference sequence with those and thus many of the minimizers will be pushed up to higher taxonomic levels. This effect will be very different between different parts of the tree of life, depending on species richness and how much of the diversity that is sequenced. We therefore come up with the idea of “classifiable regions” (CRs), i.e. regions of the reference genome that reads can be classified to using a specific kraken database. We can’t find out which regions have minimizers in the database from Kraken, and we have to find it out some other way. What we have come up with is to simulate reads from the reference genome used in the kraken database, and then classify them using the same database, then take the reads classified correctly and mapping them back to the reference genome. We can then get the number of mapped real metagenomics reads within each such CR and plot that number against the length of the CR. We also perform deduplication, to get the true biological signal and not the PCR and optical duplicates.
## Modifications
We should create a new branch that tracks the current one for any modifications we now do.
Below is a list of different areas that needs addressing. Please ingest this prompt and make a structured plan on how to tackle all areas and all ideas in this prompt.
### Estimation of classifiable regions (CRs)
Change to a sliding window slicing of the genome instead of full-on proper read simulation. This will be much faster and should still work just as fine for this use case. We can maybe use the sliding function in seqkit. We should have a discussion if it’s enough to get single-end reads for this step, or if we need 150x2 bp paired-end reads. Our metagenomics data is 150x2 bp in raw format (but can be truncated during preprocessing).
### Checkpoint
We should remove the checkpoint or drastically rework the conditions, currently it’s way to easy for a species not to be validated due to few reads / low coverages.
### Choice of mapper
Must change mapper to something that can handle extremely long contigs, the largest genomes, and one that preferably runs fast and with low memory footprint. I believe that the current mapper BBMap only allows for contigs of length 500 Mbp. I would like you to double check that.
### Reference contiguity
The pipeline must accommodate both highly contiguous genomes (chromosome level) as well as highly fragmented genomes (thousands or even millions of contigs). It will be difficult to get a good correlation for an organism with a very high-quality reference genome if it consists of only a handful chromosomes. We should discuss how we change the logic for getting the correlation. We could split the reference genomes into smaller pieces for example, but they would have to be of different size so we can get a correlation. If we split, we could also do that theoretically after mapping – i.e., no need to change the fasta file. I don’t know exactly how to approach this, we have to discuss.
### Include a read extraction pipeline and genome download
The pipeline should include extraction of reads using kraken classification output and fastq files, see “../ReadExtractorSMK”. Currently the user has to run that pipeline first, then concatenate the reads and move them into the input folder here. The user could specify a tax_id and a reference sequence accession and the pipeline then extracts the reads from the metagenomics fastq files. The pipeline can also download the genomes using ncbi datasets CLI. The user can supply an api-key to use.
### Sample-aware mapping
We want to be able to get how many (deduped) reads that map to the reference genome for each sample (week) and for each species. This way we can get new “clean” signal for each species. You can look in the folder ../six_birds to get inspiration on how to do sample-aware mapping. 
### Dynamic SLURM resource usage
For jobs that one can expect to require different amount of resources (time/memory/CPU cores), the pipeline should make an initial reasonable guess and then incrementally increase the resources across the retries (retries should be 3-5). I’m thinking of jobs that include indexing/mapping/classification for example.
We should tailor this to use on the Kebnekaise cluster of HPC2N.
### Tidy pipeline
The pipeline should make sure to mark intermediate files as temporary so that they are removed properly and don’t take up too much space.
I would also gladly appreciate suggestions on how to structure the outputs of the pipeline. Maybe we need to change it, I’m not sure.
It should be easy for a user to go into the output structure and find whatever species-related results that he or she wants, without having to look for a specific taxID.
### Output report
Final report should be a self-contained html file that is sleek, clear, modern-looking. Again, I gladly receive suggestions on how to revamp the report and how it should present the validation results.
### Taxonomy
We should use JolTax for any taxonomic work if needed. You can look in ../../../../03-DEVEL to find JolTax.
### Sanity-check
I want you to take on critical glasses and do a thorough transillumination of the whole pipeline, and really think about if we can restructure it in any way. Are all job steps necessary? How should we actually write the Snakemake code – currently we are heavily using templating. Are we using python scripts a lot? Is that the best approach? Maybe we can use more bash and pipes and python code in the rules themselves? If you were given the task to write this pipeline from scratch, how would you structure it?
### Readme
The readme should be updated and reworked to be highly user-friendly and not use too much technical jargon. It should be short and concise.

