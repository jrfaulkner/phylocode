# Working directory, speciy number of copies of analyses, initial seed
setwd("~/git_repos/jfaulkner_phylocode/Rev/")
n_replicate_analyses <- 2
first_seed <- 42

datasets <- c("bottleneck","mex_hat","broken_exponential","GP")

# For setting priors
clock_rates <- scan("tests/setup/simulated_clock_rates.txt")
const_pops <- log(scan("tests/setup/const_pop_mle.txt"))
grid_ends <- scan("tests/setup/grid_ends.txt")

# Full model (variable-tree) analyses

GMRF_template <- readLines("tests/setup/test_GMRF_bison_template.Rev")
HSRF_template <- readLines("tests/setup/test_HSRF_bison_template.Rev")

seed <- first_seed

for (d in 1:4) {
  DS <- datasets[d]
  # Generate GMRF analysis files
  for (i in 1:n_replicate_analyses) {
    this_GMRF <- GMRF_template
    this_GMRF <- gsub("THISSEED",seed,this_GMRF)
    this_GMRF <- gsub("THISREP",i,this_GMRF)
    this_GMRF <- gsub("THISALN",DS,this_GMRF)
    this_GMRF <- gsub("THISCLOCK",clock_rates[d],this_GMRF)
    this_GMRF <- gsub("THISGRID",grid_ends[d],this_GMRF)
    this_GMRF <- gsub("THISMU",const_pops[d],this_GMRF)
    cat(this_GMRF,file=paste0("tests/full_model/src/test_",DS,"_GMRF_",i,".Rev"),sep="\n")
    seed <- seed + 1
  }
  
  # Generate HSRF analysis files
  for (i in 1:n_replicate_analyses) {
    this_HSRF <- HSRF_template
    this_HSRF <- gsub("THISSEED",seed,this_HSRF)
    this_HSRF <- gsub("THISREP",i,this_HSRF)
    this_HSRF <- gsub("THISALN",DS,this_HSRF)
    this_HSRF <- gsub("THISCLOCK",clock_rates[d],this_HSRF)
    this_HSRF <- gsub("THISGRID",grid_ends[d],this_HSRF)
    this_HSRF <- gsub("THISMU",const_pops[d],this_HSRF)
    cat(this_HSRF,file=paste0("tests/full_model/src/test_",DS,"_HSRF_",i,".Rev"),sep="\n")
    seed <- seed + 1
  }
}

# Fixed tree analyses

GMRF_template <- readLines("tests/setup/test_GMRF_tree_only_template.Rev")
HSRF_template <- readLines("tests/setup/test_HSRF_tree_only_template.Rev")

for (d in 1:4) {
  DS <- datasets[d]
  # Generate GMRF analysis files
  for (i in 1:n_replicate_analyses) {
    this_GMRF <- GMRF_template
    this_GMRF <- gsub("THISSEED",seed,this_GMRF)
    this_GMRF <- gsub("THISREP",i,this_GMRF)
    this_GMRF <- gsub("THISALN",DS,this_GMRF)
    this_GMRF <- gsub("THISGRID",grid_ends[d],this_GMRF)
    this_GMRF <- gsub("THISMU",const_pops[d],this_GMRF)
    cat(this_GMRF,file=paste0("tests/tree_only/src/test_",DS,"_GMRF_",i,".Rev"),sep="\n")
    seed <- seed + 1
  }
  
  # Generate HSRF analysis files
  for (i in 1:n_replicate_analyses) {
    this_HSRF <- HSRF_template
    this_HSRF <- gsub("THISSEED",seed,this_HSRF)
    this_HSRF <- gsub("THISREP",i,this_HSRF)
    this_HSRF <- gsub("THISALN",DS,this_HSRF)
    this_HSRF <- gsub("THISGRID",grid_ends[d],this_HSRF)
    this_HSRF <- gsub("THISMU",const_pops[d],this_HSRF)
    cat(this_HSRF,file=paste0("tests/tree_only/src/test_",DS,"_HSRF_",i,".Rev"),sep="\n")
    seed <- seed + 1
  }
}