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
