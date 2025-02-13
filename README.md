<img src="doc/Logo.png" alt="decisivator logo" width=300px/> 

#### Decisivator

##### Overview
Calculate phylogenetic 'decisiveness' sensu Sanderson, Steel, and McMahon. Calculations can 
involve matrix occupancy alone (for so-called 'complete decisiveness'). The combinatorics 
involved can be prohibitive with large taxon-gene matrices. As an alternative, decisiveness 
can be approximated through the analysis of many random tree topologies. Additionally, the 
decisiveness of a user tree (or distribution of trees) can be quantified.

##### Compilation
To compile, type the following in a unix prompt:

	make

By default, this compiles using OPENMP to enable multithreaded execution for batch analyses.
If (for some reason) you wish to turn off multithreading, recompile as:

	make clean
	make mt=F

##### Usage

To run, type:

	./Decisivator [-d data_file] [-m taxon-gene_matrix] [-t tree_file] [-b burnin]
	[-n thinning] [-w taxon_weights] [-l locus_weights]

where:

	'data_file' is a simple 'vanilla' Nexus file containing sequences and defined CHARSETs.
	  - PLEASE NOTE! Only simple CHARSETs are currently supported.
	    - e.g. contiguous (X-Y) or interval (e.g. codon: X-Y\3) data are fine.
	    - CHARSET referencing is NOT allowed at present (but will be!)

	'taxon-gene_matrix' is a (legacy format) table listing taxa (rows) and genes (columns).
	  '1' indicates cell has been sequenced, while '0' indicates it has not.
	  First row should give locus names. First column should give taxon names.

	'tree_file' contains user tree(s) in Nexus format to evaluate decisiveness upon. Tree(s)
	  must be fully bifurcating.

	'burnin' is the number of trees to ignore. Only makes sense with a distribution of trees.

	'thinning' is the interval between sampling trees (i.e. where every nth tree sample will
	  be retained). Only makes sense with a distribution of trees.

	'taxon_weights' is a two-column (taxon, weight; with headers) file listing weights for
	  taxa based on some arbitrary accessibility criterion.

	'locus_weights' is a two-column (locus, weight; with headers) file for locus weights.

	NOTE: The taxon and locus weight files need not be complete. All weights are 1.0 by
	  default. Enter only weights which should be changed from the default. Or do this kind
	  of stuff within the program itself.

For help, type:

	./Decisivator -h
