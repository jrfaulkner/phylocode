# NOTE: This file is setup to be sourced from the Rev folder
# Specify number of copies of analyses (chains), initial seed
n_replicate_analyses <- 4
first_seed_1g <- 101
first_seed_1h <- 201

first_seed_2g <- 301
first_seed_2h <- 401

GMRF_1_template <- readLines("analyses/setup/hcv_GMRF_order_1_template.Rev")
HSMRF_1_template <- readLines("analyses/setup/hcv_HSMRF_order_1_template.Rev")
GMRF_2_template <- readLines("analyses/setup/hcv_GMRF_order_2_template.Rev")
HSMRF_2_template <- readLines("analyses/setup/hcv_HSMRF_order_2_template.Rev")


dir.create("analyses/hcv/src", recursive=TRUE, showWarnings=FALSE)

### ---- Order 1 -------------------------
seed_1 <- first_seed_1g
# Generate GMRF-1 analysis files
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_1_template
  this_GMRF <- gsub("THISSEED",seed_1,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  cat(this_GMRF,file=paste0("analyses/hcv/src/analyze_hcv_GMRF_1_chain_",i,".Rev"),sep="\n")
  seed_1 <- seed_1 + 1
}

# Generate HSMRF-1 analysis files
seed_1 <- first_seed_1h
for (i in 1:n_replicate_analyses) {
  this_HSMRF <- HSMRF_1_template
  this_HSMRF <- gsub("THISSEED",seed_1,this_HSMRF)
  this_HSMRF <- gsub("THISREP",i,this_HSMRF)
  cat(this_HSMRF,file=paste0("analyses/hcv/src/analyze_hcv_HSMRF_1_chain_",i,".Rev"),sep="\n")
  seed_1 <- seed_1 + 1
}

### ---- Order 2 -------------------------

# Generate GMRF-2 analysis files
seed_2 <- first_seed_2g
for (i in 1:n_replicate_analyses) {
  this_GMRF <- GMRF_2_template
  this_GMRF <- gsub("THISSEED",seed_2,this_GMRF)
  this_GMRF <- gsub("THISREP",i,this_GMRF)
  cat(this_GMRF,file=paste0("analyses/hcv/src/analyze_hcv_GMRF_2_chain_",i,".Rev"),sep="\n")
  seed_2 <- seed_2 + 1
}

# Generate HSMRF-2 analysis files
seed_2 <- first_seed_2h
for (i in 1:n_replicate_analyses) {
  this_HSMRF <- HSMRF_2_template
  this_HSMRF <- gsub("THISSEED",seed_2,this_HSMRF)
  this_HSMRF <- gsub("THISREP",i,this_HSMRF)
  cat(this_HSMRF,file=paste0("analyses/hcv/src/analyze_hcv_HSMRF_2_chain_",i,".Rev"),sep="\n")
  seed_2 <- seed_2 + 1
}


