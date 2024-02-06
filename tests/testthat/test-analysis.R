library(tidyverse)
set.seed(123)

# ------------------------------------------------------------------------------
# helper functions
# ------------------------------------------------------------------------------
computeCounts = function(population, condition, mean, sd, n, parent_count = NULL) {

  pct_parent = rnorm(n, mean = mean, sd = sd)

  if (is.null(parent_count)) {
    parent_count = 10000
  }

  return(data.frame(name = paste0(condition, "_", 1:n),
                    Population = population,
                    condition = condition,
                    Count = round(parent_count*pct_parent),
                    ParentCount = parent_count,
                    PercentParent = pct_parent
  ))
}


runModifierAnalysis = function(counts) {
  summary = summarizePopulations(counts, verbose = F)
  summary = processGatingPaths(summary, max_depth = 2)


  stats = generateStatistics(counts, summary, "healthy", "disease", verbose = F)
  stats = processGatingPaths(stats, max_depth = 2)

  report = modifierAnalysis(stats, counts, p_cutoff = 0.05, verbose = F)

  return(report$path)
}

# ------------------------------------------------------------------------------
# Setup parameters that are used across all tests
# ------------------------------------------------------------------------------
n = 30 # number of samples that

# ------------------------------------------------------------------------------
# 1. test: only parent population is significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+"

expect = c("/R+/A+")

counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n)
)


test_that("modifier analysis: test 1", {
  expect_equal(runModifierAnalysis(counts), expect)
})

# ------------------------------------------------------------------------------
# 2. test: parent populations are not significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+/B+"
expect = c("/R+/A+/B+")

counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n),
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n)
)


test_that("modifier analysis: test 2", {
  expect_equal(runModifierAnalysis(counts), expect)
})

# ------------------------------------------------------------------------------
# 3. test: one parent populations has the same distribution as the modifier
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/B+"
expect = c("/R+/B+")

counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.6, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)


test_that("modifier analysis: test 3", {
  expect_equal(runModifierAnalysis(counts), expect)
})

# ------------------------------------------------------------------------------
# 4. test: all populations are different and significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+, "/R+/B+, "/R+/A+/B+"
expect = sort(c("/R+/A+/B+", "/R+/B+", "/R+/A+"))

counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.1, 0.1, n), # same as modifier
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.6, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)

test_that("modifier analysis: test 4", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 5. test: one parent and modifier populations are different and significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+, "/R+/A+/B+"
expect = sort(c("/R+/A+", "/R+/A+/B+"))

counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.1, 0.1, n), # same as modifier
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.2, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)

test_that("modifier analysis: test 5", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})


# ------------------------------------------------------------------------------
# Update test function to run to depth 3
# ------------------------------------------------------------------------------
runModifierAnalysis = function(counts) {
  summary = summarizePopulations(counts, verbose = F)
  summary = processGatingPaths(summary, max_depth = 3)


  stats = generateStatistics(counts, summary, "healthy", "disease", verbose = F)
  stats = processGatingPaths(stats, max_depth = 3)

  report = modifierAnalysis(stats, counts, p_cutoff = 0.05, verbose = F)

  return(report$path)
}

# ------------------------------------------------------------------------------
# 6. test: three gates, only one is significant - variation 1
# ------------------------------------------------------------------------------
expect = sort(c("/R+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 6", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 7. test: three gates, only one is significant - variation 2
# ------------------------------------------------------------------------------
expect = sort(c("/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 7", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 8. test: three gates, only one is significant - variation 3
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 8", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})


# ------------------------------------------------------------------------------
# 9. test: parent and leafs are significant, but not significant vs each other
# ------------------------------------------------------------------------------
expect = sort(c("/R+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 9", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 10. test: parent and leafs are significant, and significant across each other
# ------------------------------------------------------------------------------
expect = sort(c("/R+/C+", "/R+/A+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.3, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 10", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})


# ------------------------------------------------------------------------------
# 11. test: parent is not significant, leafs have same distribution
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 11", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})


# ------------------------------------------------------------------------------
# 12. test: parent is not significant, leafs have different distribution
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/C+", "/R+/A+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.5, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 12", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 13. test: all are significant, all are different
# ------------------------------------------------------------------------------
expect = sort(c("/R+/C+", "/R+/A+/C+", "/R+/B+/C+", "/R+/A+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.3, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.1, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.5, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.1, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.6, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 13", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 14. test: only level 1 are significant
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/C+", "/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 14", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})


# ------------------------------------------------------------------------------
# 15. test: level 2 are significant but different, and level 3 is significant too
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/B+/C+", "/R+/A+/C+", "/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 15", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 16. test: level 2 are significant but different, and level 3 is significant too, all different logFCs
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/B+/C+", "/R+/A+/C+", "/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.5, 0.05, n)
)

test_that("modifier analysis: test 16", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 17. test:
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/B+/C+", "/R+/A+/C+", "/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.2, 0.05, n)
)

test_that("modifier analysis: test 17", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})

# ------------------------------------------------------------------------------
# 18. test:
# ------------------------------------------------------------------------------
expect = sort(c("/R+/A+/B+/C+", "/R+/A+/C+", "/R+/B+/C+"))

counts = rbind(
  computeCounts("/R+/C+", "healthy", 0.4, 0.05, n),
  computeCounts("/R+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/C+", "healthy", 0.8, 0.05, n),
  computeCounts("/R+/A+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "healthy", 0.2, 0.05, n),
  computeCounts("/R+/B+/C+", "disease", 0.2, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "healthy", 0.5, 0.05, n),
  computeCounts("/R+/A+/B+/C+", "disease", 0.1, 0.05, n)
)

test_that("modifier analysis: test 18", {
  expect_equal(sort(runModifierAnalysis(counts)), expect)
})
