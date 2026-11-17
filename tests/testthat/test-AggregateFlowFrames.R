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
                 'file2.fcs: FL2-A is NA ("CD4" in first file)',
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

# ---- Files with different channels ----
# A flowSet cannot hold files with different channels, so these tests write
# small FCS files. Every channel X gets the marker "M_X".
write_files <- function(...){
  dir <- tempfile("agg")
  dir.create(dir)
  channelSets <- list(...)
  vapply(seq_along(channelSets), function(i){
    channels <- channelSets[[i]]
    m <- matrix(stats::runif(50 * length(channels), 1, 5),
                ncol = length(channels), dimnames = list(NULL, channels))
    ff <- flowCore::flowFrame(m)
    p <- flowCore::parameters(ff)
    p$desc <- paste0("M_", channels)
    flowCore::parameters(ff) <- p
    path <- file.path(dir, paste0("file", i, ".fcs"))
    flowCore::write.FCS(ff, path)
    path
  }, character(1))
}

test_that("missing channels are filled with zeros by default, with details", {
  files <- write_files(c("A", "B", "C"), c("A", "B"))
  w <- testthat::capture_warnings(
    ff <- AggregateFlowFrames(files, 40, silent = TRUE))
  expect_length(w, 1) # no extra marker warning for the zero-filled channel
  expect_match(w, "file2.fcs: C", fixed = TRUE)
  expect_match(w, "filled with zeros", fixed = TRUE)
  e <- flowCore::exprs(ff)
  expect_true(all(e[e[, "File"] == 2, "C"] == 0))
  expect_true(all(e[e[, "File"] == 1, "C"] != 0))
})

test_that("missingChannels = 'drop' keeps only shared channels", {
  files <- write_files(c("A", "B", "C"), c("A", "B"))
  expect_warning(
    ff <- AggregateFlowFrames(files, 40, silent = TRUE,
                              missingChannels = "drop"),
    "left out of the aggregate: C", fixed = TRUE)
  expect_false("C" %in% colnames(ff))
  expect_true(all(c("A", "B") %in% colnames(ff)))
})

test_that("missingChannels = 'error' stops and lists the missing channels", {
  files <- write_files(c("A", "B", "C"), c("A", "B"), c("A", "C"))
  expect_error(AggregateFlowFrames(files, 40, silent = TRUE,
                                   missingChannels = "error"),
               "file2.fcs: C\n  file3.fcs: B", fixed = TRUE)
})

test_that("channels missing from the first file are reported", {
  files <- write_files(c("A", "B"), c("A", "B", "C"))
  expect_warning(ff <- AggregateFlowFrames(files, 40, silent = TRUE),
                 "not in the first file (file1.fcs)", fixed = TRUE)
  expect_false("C" %in% colnames(ff))
})

test_that("channels given as marker names do not give a false warning", {
  files <- write_files(c("A", "B", "C"), c("A", "B", "C"))
  expect_no_warning(
    ff <- AggregateFlowFrames(files, 40, silent = TRUE,
                              channels = c("M_A", "M_C")))
  expect_equal(setdiff(colnames(ff), c("File", "File_scattered",
                                       "Original_ID")),
               c("A", "C"))
})

test_that("requested channels missing in later files are reported", {
  files <- write_files(c("A", "B", "C"), c("A", "B"))
  expect_warning(AggregateFlowFrames(files, 40, silent = TRUE,
                                     channels = c("M_A", "M_C")),
                 "the requested channels:\n  file2.fcs: C", fixed = TRUE)
})
