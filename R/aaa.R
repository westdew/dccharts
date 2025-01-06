the <- new.env(parent = emptyenv())
the$cache <- list()

cache_set <- function(key, val) {
  the$cache[[key]] <- val
}

cache_get <- function(key) {
  the$cache[[key]]
}
