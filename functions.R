summarizePopulations = function(counts) {

  message(paste0("Summarizing ", length(unique(counts$Population)), " unique populations."))

  summary = counts %>%
    dplyr::select(name, path = Population, condition, Count, ParentCount, PercentParent)%>%
    dplyr::filter(!is.na(PercentParent)) %>%
    group_by(path, condition) %>%
    summarise(mean = mean(PercentParent),
              sd = sd(PercentParent),
              median = median(PercentParent),
              nCells = sum(Count))

  # we need to remove the grouping, or it will later slow down other analysis steps
  return(as.data.frame(summary))
}

processGatingPaths = function(summary, max_depth = 1) {

  columns = colnames(summary)

  if (!"path" %in% columns) {
    summary$path = summary$Population
  }

  # determine gating depth for each population
  summary$depth = str_count(summary$path, "/")

  # determine gating depth of the root path
  root_path_length = max(summary$depth)-max_depth

  # split all gating paths into a matrix so that gating results can be more
  # easily extracted later
  pathing = str_split(str_remove(summary$path, '^/'), "/", simplify = TRUE)

  # determine the root path
  root_path = paste(pathing[nrow(summary), 1:root_path_length], collapse = "/")
  summary$root_node = root_path

  # determine terminal node
  summary$terminal_node = pathing[cbind(1:nrow(pathing), summary$depth)]

  # determine gating depth after root path
  summary$depth = summary$depth - root_path_length

  # No need to extract modifiers if max_depth is 1
  if (max_depth == 1) {return(summary)}

  # determine what the maximum number of modifiers must have been
  max_modifiers = max_depth -1

  # determine all modifiers and add columns for them, will generate columns of
  # the type modifier_1, modifier_2, ... up to the maximum number of modifiers
  for (i_modifier in 1:max_modifiers ){
    column_name = paste0("modifier", "_", i_modifier)

    summary[[column_name]] = ""

    idx = which(summary$depth>i_modifier)

    summary[idx, column_name] = pathing[cbind(idx, root_path_length+i_modifier)]
  }

  return(summary)
}

generateStatistics = function(counts, summary, control, studyGroups) {
  paths = unique(summary$path)
  nPops = length(paths)
  combined.stats <- vector("list", length = nPops)

  message(paste0("Calculating statistics for ", nPops, " unique populations"))
  i = 0
  pb = txtProgressBar(min = i, max = nPops, style = 3)

  # calculate statistics by looping over populations
  # This could probably get turned into a summary/group_by step?
  for (ipath in 1:nPops) {
    i = i+1
    setTxtProgressBar(pb, i)

    path_ = paths[i]
    # print(path_)
    stats = c()

    # get the pct parents for that population
    x = dplyr::filter(counts, condition == control, Population == path_)$PercentParent
    x = x[!is.na(x)] # but filter out the entries that are NA - this is the case when the parent already had 0 cells for that sample

    x_counts = dplyr::filter(counts, condition == control, Population == path_)$Count

    for (group_ in studyGroups ) {
      y = dplyr::filter(counts, condition == group_, Population == path_)$PercentParent
      y = y[!is.na(y)]

      y_counts = dplyr::filter(counts, condition == group_, Population == path_)$Count

      if (length(x) == 0 | length(y) == 0) {
        # if one of the study groups has just NA's for the pct parent counts
        log2FC = NA
        p.value = NA
      } else {
        log2FC = log2(median(y)) - log2(median(x))
        test = wilcox.test(x, y, exact = FALSE)
        p.value = test$p.value
      }

      stats_ = data.frame(path = path_,
                          condition = group_,
                          log2Fold = log2FC,
                          pvalue = p.value,
                          group_1_name = control,
                          group_2_name = group_,
                          group_1_median_pct_parent = median(x),
                          group_2_median_pct_parent = median(y),
                          group_1_sd_pct_parent = sd(x),
                          group_2_sd_pct_parent = sd(y),
                          group_1_median_count = median(x_counts),
                          group_2_median_count = median(y_counts),
                          group_1_sd_count = sd(x_counts),
                          group_2_sd_count = sd(y_counts),
                          group_1_total_cells = sum(x_counts),
                          group_2_total_cells = sum(y_counts))

      stats = rbind(stats, stats_)
    }

    combined.stats[[ipath]] = stats
  }

  close(pb)

  combined.stats.df = do.call(rbind, combined.stats)

  combined.stats.df = dplyr::filter(combined.stats.df, !is.na(pvalue))

  return(combined.stats.df)

}

modifierAnalysis = function(fullstats, counts, p_cutoff = 0.05) {

  group_comparisons = distinct(dplyr::select(fullstats, group_1_name, group_2_name))

  full_report = c()
  removed_from_report = c()

  for (icomparison in 1:nrow(group_comparisons)) {
    control_group = group_comparisons$group_1_name[icomparison]
    study_group = group_comparisons$group_2_name[icomparison]

    message(paste0("Starting modifier analysis of ", control_group, " vs ", study_group))

    stats = dplyr::filter(fullstats, condition == study_group)
    report = c()

    max_depth = max(stats$depth)

    # Find all the nodes that we have
    # First remove + and -, find uniques, and then readd +/- in regexp compatible form
    nodes = unique(str_remove(dplyr::filter(stats, depth > 0)$terminal_node, "\\+|\\-"))
    nodes = c(paste0(nodes, "\\+"), paste0(nodes, "\\-"))

    # Extract rood node, and turn it into a regexp compatible version
    root_node = stats$root_node[1]
    root_node = str_replace_all(root_node, "\\+", "\\\\+")
    root_node = str_replace_all(root_node, "\\-", "\\\\-")
    root_node = paste0("/", root_node)

    # Add a column modifier_path by removing the root node and also removing the terminal node
    stats$modifier_path = sub(root_node, "", stats$path)
    stats$modifier_path = sub("/[^/]*$", "", stats$modifier_path)

    # Build a binary mask for all modifier paths
    modifier_masks = lapply(stats$modifier_path, function(x) { ifelse(str_detect(x, nodes), "1", "0") })
    for (i in 1:length(modifier_masks)) {
      stats$modifier_mask[i] = unlist(paste0(modifier_masks[[i]], collapse = ""))
    }

    subpopulations = list()
    for (depth_ in max_depth:1) {
      subpopulations[[depth_]] = dplyr::filter(stats, depth == depth_, pvalue < p_cutoff)
    }

    for (depth_ in max_depth:1) {
      nPops = nrow(subpopulations[[depth_]])

      message(paste0("Modifier analysis of ", nPops, " subpopulations with p-value < ", p_cutoff, " at gating depth ", depth_, "."))
      iPop = 0

      if (nPops == 0) {next}

      pb = txtProgressBar(min = iPop, max = nPops, style = 3)

      for (iPop in 1:nPops) {
        setTxtProgressBar(pb, iPop)

        population_ = subpopulations[[depth_]][iPop, ]

        # Turn all the modifiers into a vector
        modifier_vec = unlist(strsplit(str_remove(population_$modifier_path, '^/'), "/"))

        if (depth_ == 1) {
          report = rbind(report, population_)
          next
        }

        nmodifier_ = depth_-1

        # Calculate the combinations of the modifiers
        if (nmodifier_ == 1) {
          combinations = matrix("", ncol = 1)
        } else {
          combinations = combn(modifier_vec, m = nmodifier_-1)
        }

        background_corrected_x = dplyr::filter(counts, condition == study_group, Population == population_$path, !is.na(PercentParent))$PercentParent -
          median(dplyr::filter(counts, condition == control_group, Population == population_$path, !is.na(PercentParent) )$PercentParent)

        case_counter = matrix(0,
                              nrow = 2, ncol = 2,
                              dimnames = list(c("non_significant_FDR_in_parent", "significant_FDR_in_parent"),
                                              c("no_significance_between_node_and_parent", "significance_between_node_and_parent")
                              ))
        parents = c()
        for (iCombination in 1:ncol(combinations)) {
          combination_mask = paste0(ifelse(str_detect(paste0(combinations[ , iCombination], collapse = "."), nodes), "1", "0"), collapse = "")
          parent_ = dplyr::filter(stats, modifier_mask == combination_mask, terminal_node == population_$terminal_node, condition == study_group)

          if (nrow(parent_) == 0) {next}

          parents = rbind(parents, parent_)

          background_corrected_y = dplyr::filter(counts, condition == study_group, Population == parent_$path, !is.na(PercentParent))$PercentParent -
            median(dplyr::filter(counts, condition == control_group, Population == parent_$path, !is.na(PercentParent) )$PercentParent)

          p_value_node_vs_parent = wilcox.test(background_corrected_x, background_corrected_y, exact = FALSE)$p.value

          # if R/A/B/C is investigated
          # if there is a significant difference between R/A/C and R/A/B/C
          if (p_value_node_vs_parent < p_cutoff) {
            # R/A/C vs R/A/B/C IS significant
            col = "significance_between_node_and_parent"
          } else {
            col = "no_significance_between_node_and_parent"
          }

          if (parent_$pvalue < p_cutoff) {
            # If R/A/C by itself is significant
            row = "significant_FDR_in_parent"
          } else {
            row = "non_significant_FDR_in_parent"
          }

          case_counter[row, col] =+ 1

        }

        if (sum(case_counter["significant_FDR_in_parent", ]) == 1 &
            case_counter["significant_FDR_in_parent", "no_significance_between_node_and_parent"] == 1 ) {
          # This is the case where one parent has a significant FDR, and the select subpopulation compared to that parent has no significant difference
          # In this case, we can assume that the observed difference in the node is driven by the parent.
          # We need to make sure that a) there is ONLY ONE parent that is significant and b) that this parent shows no difference in the distribution of values between itself and the subpopulation
          # We only want to report the parent, and remove the node from the report (in case it has been added previously by a higher gating depth)

          if (!is.null(report)) {
            report = dplyr::filter(report, path != population_$path)
          }
          # The parent is getting reported outside the if condition

          # For debugging reasons, I am also keeping track of subpopulations that get removed
          removed_from_report = rbind(removed_from_report, population_)
        } else {
          report = rbind(report, population_)
        }

        # report all parents that have a significant FDR (if there are parents at all)
        if (!is.null(parents)) {
          report = rbind(report, dplyr::filter(parents, pvalue < p_cutoff) )
        }

      }
      close(pb)
    }

    report = dplyr::select(report, -modifier_mask, -modifier_path)
    report = distinct(report)

    full_report = rbind(full_report, report)

    if (!is.null(removed_from_report)) {
      removed_from_report = distinct(removed_from_report)
    }

  }
  return(full_report)
}
