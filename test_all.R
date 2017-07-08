library(methods)
library(devtools)
library(testthat)

if (interactive()) {
  load_all(".")
} else {
  library(mcMST)
}

test_dir("tests/testthat")
