# Small helper: a flowFrame with 3 channels and the given marker names
# (NA = channel without a marker).
make_ff <- function(markers, n = 50){
  channels <- c("FSC-A", "FL1-A", "FL2-A")
  m <- matrix(stats::rnorm(n * 3), ncol = 3,
              dimnames = list(NULL, channels))
  ff <- flowCore::flowFrame(m)
  p <- flowCore::parameters(ff)
  p$desc <- c(NA, markers)
  flowCore::parameters(ff) <- p
  ff
}

aggregate_set <- function(...){
  fs <- flowCore::flowSet(list(...))
  flowCore::sampleNames(fs) <- paste0("file", seq_along(list(...)), ".fcs")
  AggregateFlowFrames(fs, cTotal = 40, silent = TRUE)
}

test_that("example from the documentation still works", {
  fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
  expect_no_warning(
    ff <- AggregateFlowFrames(c(fileName, fileName), 1000, silent = TRUE))
  expect_equal(unname(nrow(ff)), 1000)
  expect_true(all(c("File", "File_scattered", "Original_ID") %in%
                    colnames(ff)))
})

test_that("no warning when all files have the same markers", {
  expect_no_warning(aggregate_set(make_ff(c("CD3", "CD4")),
                                  make_ff(c("CD3", "CD4"))))
})

test_that("renamed marker is named in the warning", {
  expect_warning(aggregate_set(make_ff(c("CD3", "CD4")),
                               make_ff(c("CD3e", "CD4"))),
                 'file2.fcs: FL1-A is "CD3e" ("CD3" in first file)',
                 fixed = TRUE)
})

test_that("missing marker is named in the warning", {
  expect_warning(aggregate_set(make_ff(c("CD3", "CD4")),
                               make_ff(c("CD3", NA))),
                 'file2.fcs: FL2-A is <none> ("CD4" in first file)',
                 fixed = TRUE)
})

test_that("swapped markers are detected", {
  expect_warning(aggregate_set(make_ff(c("CD3", "CD4")),
                               make_ff(c("CD4", "CD3"))),
                 "FL1-A is \"CD4\".*FL2-A is \"CD3\"")
})

test_that("mismatches from several files end up in one warning", {
  w <- testthat::capture_warnings(
    aggregate_set(make_ff(c("CD3", "CD4")),
                  make_ff(c("CD3e", "CD4")),
                  make_ff(c("CD3", "CD8"))))
  expect_length(w, 1)
  expect_match(w, "file2.fcs", fixed = TRUE)
  expect_match(w, "file3.fcs", fixed = TRUE)
})
