# Working directory, speciy number of copies of analyses, initial seed
setwd("~/git_repos/jfaulkner_phylocode/Rev/")
n_replicate_analyses <- 2
first_seed <- 42

# Variable-tree analyses

GMRF_template <- readLines("tests/setup/test_GMRF_bison_template.Rev")
HSRF_template <- readLines("tests/setup/test_HSRF_bison_template.Rev")

seed <- first_seed

# constant-pop dataset
dataset <- "constant_population"

# Generate GMRF analysis files
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_template
  this_GMRF <- gsub("THISSEED",seed,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  this_GMRF <- gsub("THISALN",dataset,this_GMRF)
  cat(this_GMRF,file=paste0("tests/src/test_constant_population_GMRF_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# Generate HSRF analysis files
for (i in 1:n_replicate_analyses) {
  this_HSRF <- HSRF_template
  this_HSRF <- gsub("THISSEED",seed,this_HSRF)
  this_HSRF <- gsub("THISREP",i,this_HSRF)
  this_HSRF <- gsub("THISALN",dataset,this_HSRF)
  cat(this_HSRF,file=paste0("tests/src/test_constant_population_HSRF_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# epoch dataset
dataset <- "epoch"

# Generate GMRF analysis files
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_template
  this_GMRF <- gsub("THISSEED",seed,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  this_GMRF <- gsub("THISALN",dataset,this_GMRF)
  cat(this_GMRF,file=paste0("tests/src/test_epoch_GMRF_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# Generate HSRF analysis files
for (i in 1:n_replicate_analyses) {
  this_HSRF <- HSRF_template
  this_HSRF <- gsub("THISSEED",seed,this_HSRF)
  this_HSRF <- gsub("THISREP",i,this_HSRF)
  this_HSRF <- gsub("THISALN",dataset,this_HSRF)
  cat(this_HSRF,file=paste0("tests/src/test_epoch_HSRF_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# Fixed tree analyses

GMRF_template <- readLines("tests/setup/test_GMRF_tree_only_template.Rev")
HSRF_template <- readLines("tests/setup/test_HSRF_tree_only_template.Rev")

dataset <- "constant_population"

# Generate GMRF analysis files
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_template
  this_GMRF <- gsub("THISSEED",seed,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  this_GMRF <- gsub("THISALN",dataset,this_GMRF)
  cat(this_GMRF,file=paste0("tests/src/test_constant_population_GMRF_tree_only_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# Generate HSRF analysis files
for (i in 1:n_replicate_analyses) {
  this_HSRF <- HSRF_template
  this_HSRF <- gsub("THISSEED",seed,this_HSRF)
  this_HSRF <- gsub("THISREP",i,this_HSRF)
  this_HSRF <- gsub("THISALN",dataset,this_HSRF)
  cat(this_HSRF,file=paste0("tests/src/test_constant_population_HSRF_tree_only_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# epoch dataset
dataset <- "epoch"

# Generate GMRF analysis files
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_template
  this_GMRF <- gsub("THISSEED",seed,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  this_GMRF <- gsub("THISALN",dataset,this_GMRF)
  cat(this_GMRF,file=paste0("tests/src/test_epoch_GMRF_tree_only_",i,".Rev"),sep="\n")
  seed <- seed + 1
}

# Generate HSRF analysis files
for (i in 1:n_replicate_analyses) {
  this_HSRF <- HSRF_template
  this_HSRF <- gsub("THISSEED",seed,this_HSRF)
  this_HSRF <- gsub("THISREP",i,this_HSRF)
  this_HSRF <- gsub("THISALN",dataset,this_HSRF)
  cat(this_HSRF,file=paste0("tests/src/test_epoch_HSRF_tree_only_",i,".Rev"),sep="\n")
  seed <- seed + 1
}


