## List of generics
##- DIST
##- LOG
##- EXP
##- PARALLEL
##- UNROLL
##- UNWRAP
##- FRECHET
##- PGA

DIST <- function(a,b){
  UseMethod("DIST")
}

LOG <- function(p,v){
  UseMethod("LOG")
}

EXP <- function(p,v){
  UseMethod("EXP")
}

PARALLEL <- function(start, end, vec){
  UseMethod("PARALLEL")
}

UNROLL <- function(data,...){
  UseMethod("UNROLL")
}

UNWRAP <- function(data, timestamp, base){
  UseMethod("UNWRAP")
}

# WRAP <- function(data, timestamp, base){
#   UseMethod("WRAP")
# }

FRECHET <- function(data, alpha){
  UseMethod("FRECHET")
}

PGA <- function(data,...){
  UseMethod("PGA")
}

PPC <- function(data, lambda, ...){
  UseMethod("PPC")
}
