# Test modifier analysis


library(tidyverse)
# Create a hypothetical dataset by sampling from a Gaussian distribution for the
# percent parent

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

# ------------------------------------------------------------------------------
# 0. case: only parent population is significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+"
n = 30
counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n)
)

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
summary = summarizePopulations(counts)
summary = processGatingPaths(summary, max_depth = 2)


stats = generateStatistics(counts, summary, "healthy", "disease")
stats = processGatingPaths(stats, max_depth = 2)

report = modifierAnalysis(stats, counts, p_cutoff = 0.05)

print(report$path)


# ------------------------------------------------------------------------------
# 1. case: parent populations are not significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+/B+"
n = 30
counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n),
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n)
)

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
summary = summarizePopulations(counts)
summary = processGatingPaths(summary, max_depth = 2)


stats = generateStatistics(counts, summary, "healthy", "disease")
stats = processGatingPaths(stats, max_depth = 2)

report = modifierAnalysis(stats, counts, p_cutoff = 0.05)

print(report$path)

# ------------------------------------------------------------------------------
# 2. case: one parent populations has the same distribution as the modifier
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/B+"
n = 30
counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/A+", "disease", 0.4, 0.1, n),
  computeCounts("/R+/B+", "healthy", 0.6, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
summary = summarizePopulations(counts)
summary = processGatingPaths(summary, max_depth = 2)


stats = generateStatistics(counts, summary, "healthy", "disease")
stats = processGatingPaths(stats, max_depth = 2)

report = modifierAnalysis(stats, counts, p_cutoff = 0.05)

print(report$path)


# ------------------------------------------------------------------------------
# 3. case: all populations are different and significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+, "/R+/B+, "/R+/A+/B+"
n = 30
counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.1, 0.1, n), # same as modifier
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "healthy", 0.4, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.6, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
summary = summarizePopulations(counts)
summary = processGatingPaths(summary, max_depth = 2)


stats = generateStatistics(counts, summary, "healthy", "disease")
stats = processGatingPaths(stats, max_depth = 2)

report = modifierAnalysis(stats, counts, p_cutoff = 0.05)

print(report$path)

# ------------------------------------------------------------------------------
# 4. case: one parent and modifier populations are different and significant
# ------------------------------------------------------------------------------
# Expected output in report: "/R+/A+, "/R+/A+/B+"
n = 30
counts = rbind(
  computeCounts("/R+/A+", "healthy", 0.1, 0.1, n), # same as modifier
  computeCounts("/R+/A+", "disease", 0.2, 0.1, n), # same as modifier
  computeCounts("/R+/B+", "healthy", 0.2, 0.1, n),
  computeCounts("/R+/B+", "disease", 0.2, 0.1, n),
  computeCounts("/R+/A+/B+", "healthy", 0.6, 0.1, n), # same as parent
  computeCounts("/R+/A+/B+", "disease", 0.2, 0.1, n) # same as parent
)

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
summary = summarizePopulations(counts)
summary = processGatingPaths(summary, max_depth = 2)


stats = generateStatistics(counts, summary, "healthy", "disease")
stats = processGatingPaths(stats, max_depth = 2)

report = modifierAnalysis(stats, counts, p_cutoff = 0.05)

print(report$path)
