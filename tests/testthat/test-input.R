test_that("simple_xlsx", {
  expect_no_error(read_xlsx(system.file("extdata", "simple_geomx.xlsx", package = "stgeomx")))
})

test_that("duplicate_probe_xlsx", {
  expect_error(read_xlsx(system.file("extdata", "test_input_duplicate_probe.xlsx", package = "stgeomx")), "duplicate probe found")
})

test_that("simple_tsv", {
  expect_no_error(read_tsv(system.file("extdata", "simple_geomx_sample.tsv", package = "stgeomx"), system.file("extdata", "simple_geomx_count.tsv", package = "stgeomx")))
})

test_that("duplicate_probe_tsv", {
  expect_error(read_tsv(system.file("extdata", "test_input_duplicate_probe_sample.tsv", package = "stgeomx"), system.file("extdata", "test_input_duplicate_probe_count.tsv", package = "stgeomx"), "duplicate probe found"))
})

test_that("blank_count_xlsx", {
  blank_count <- read_xlsx(system.file("extdata", "blank_count.xlsx", package = "stgeomx"))
  expect_equal(blank_count[["counts"]][["s01_r001_glucagon+"]][1], 0)
})

test_that("blank_count_xlsx", {
  blank_count <- read_tsv(system.file("extdata", "blank_count_sample.tsv", package = "stgeomx"), system.file("extdata", "blank_count_count.tsv", package = "stgeomx"))
  expect_equal(blank_count[["counts"]][["s01_r001_glucagon+"]][1], 0)
})
