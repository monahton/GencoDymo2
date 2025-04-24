# Inside the file tests/testthat/test-get_latest_release.R
test_that("get_latest_release works for human species", {
  result <- get_latest_release("human", verbose = FALSE)
  expect_match(result, "^release_\\d+$")  # Check if the result matches a release format like "release_42"
})

test_that("get_latest_release throws an error for invalid species", {
  expect_error(get_latest_release("dog"), "Invalid species. Please use 'human' or 'mouse'.")
})

