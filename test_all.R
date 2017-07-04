library(methods)
library(devtools)
library(testthat)

if (interactive()) {
  load_all(".")
} else {
  library(rmoco)
}

test_dir("tests/testthat")
