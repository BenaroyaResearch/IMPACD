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


n = 30

# ------------------------------------------------------------------------------
# 0. case: only parent population is significant
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


test_that("one parent population", {
  expect_equal(runModifierAnalysis(counts), expect)
})
