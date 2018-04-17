This folder contains everything needed to run RevBayes analyses of the 2013 bison data and to perform some basic testing of the Rev implementation of the HSRF coalescent

	For stability, these can all be run with the development branch of the fork of Rev available at: https://github.com/afmagee/revbayes/tree/development

	To run any tests or analyses, Rev should be compiled from this source, and called from the top level of /analyses or /tests as appropriate.

/analyses contains the requisite scripts and files for analyzing the data

	/analyses/data contains the Bison data
	
	/analyses/output is where the output of the analyses go
	
	/analyses/setup
		
		To generate the scripts for analysis, use /analyses/setup/generate_analysis_scripts.R 
		
			The number of replicate analyses can be changed by altering this file
		
		To change any details of the HSRF or GMRF analyses, alter the templates and re-generate the analysis files

	/analyses/src contains the scripts to run the analyses
	
	
/tests contains scripts to simulate and perform basic sanity checks of the Rev implementation of the HSRF and GMRF models
	
	/setup contains scripts and datafiles to simulate test datasets and write test analysis scripts
		
		/simtraj.r is used for simulating coalescent models (need phylofuncs.r, package phylodyn)
		
		/subsample_trees_for_full_Rev_test.R makes smaller versions of the full trees to use in testing the joint/full model in Rev
		
		/sim_seqs.Rev simulates JC69 alignments on the subsampled trees (using clock rates that correspond roughly to the substitution (MrBayes-style) tree length of the Bison data
		
		/generate_test_scripts.R generates the Rev scripts for analysis, uses the misc. text files to set grid ends and priors where needed in the Rev scripts
	
	/tree_only is for testing the Rev HSRF and GMRF implementations on fixed trees
		
		/data contains simulated data	
		
		/output holds the output of the analyses
		
		/src contains the generated analysis scripts

	/full_model is for testing the Rev HSRF and GMRF implementations while inferring the tree
		
		/data contains simulated data	
		
		/output holds the output of the analyses
		
		/src contains the generated analysis scripts

