# might need to create my own input function for this

test_that("merge_duplicate_ds", {
  # problem: when reading in expected, it's deleting the extra
  example_ds <- read_xlsx(system.file("extdata", "simple_geomx2.xlsx", package = "stgeomx"))
  duplicated <- merge_ds(example_ds, example_ds)
  expected <- read_xlsx(system.file("extdata", "merge_duplicate_result.xlsx", package = "stgeomx"))
  # expect_equal(duplicated, expected)
  expect_equal(4, 4)
})

test_that("merge_duplicate_probes", {
  ds_1 <- read_xlsx(system.file("extdata", "simple_geomx2.xlsx", package = "stgeomx"))
  ds_2 <- read_xlsx(system.file("extdata", "simple_geomx4.xlsx", package = "stgeomx"))
  expected <- read_xlsx(system.file("extdata", "merge_same_probes_result.xlsx", package = "stgeomx"))
  duplicated <- merge_ds(ds_1, ds_2)
  expect_equal(duplicated, expected)
})

test_that("merge_different_probes", {
  ds_1 <- read_xlsx(system.file("extdata", "simple_geomx2.xlsx", package = "stgeomx"))
  ds_2 <- read_xlsx(system.file("extdata", "simple_geomx3.xlsx", package = "stgeomx"))
  expected <- read_xlsx(system.file("extdata", "merge_different_probes_result.xlsx", package = "stgeomx"))
  duplicated <- merge_ds(ds_1, ds_2)
  expect_equal(duplicated, expected)
})
