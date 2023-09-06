library(tidyr)
library(tidyverse)
library(readxl)

devtools::load_all()

IMPACD_folder = "/nfs/spribitzer/projects/IMPACD"
IMPACD_R_folder = "/nfs/spribitzer/projects/IMPACD/R"

# Make sure that the .RProj file that your Rstudio is using is located in the package directory.
# Else you will get some annoying "can't do this" messages from `usethis`

# ------------------------------------------------------------------------------
# Get transformers from the workspace file
# ------------------------------------------------------------------------------

# Extracting data
ws <- open_flowjo_xml(file.path(IMPACD_folder,"Original data/Surface CD4+ T cell analysis 2021-05-06.wsp"))
gs_flowJo = flowjo_to_gatingset(ws, name = "batch#1")

trans <- gh_get_transformations(gs_flowJo[[1]])
transList_wsp <- transformList(from = names(trans), trans)

# Save data
usethis::use_data(transList_wsp)

# ------------------------------------------------------------------------------
# Locate fcs files and combine with annotation
# ------------------------------------------------------------------------------
fcsLocation = file.path(IMPACD_folder, "Original data/Export CD4+ T cells")
fcsFiles = list.files(fcsLocation, full.names = T, pattern = "*.fcs")
fcsFiles = data.frame(fileName = fcsFiles)

for (iRow in 1:nrow(fcsFiles)) {
  str_ = fcsFiles$fileName[iRow]
  fcsFiles$sampleName[iRow] = str_split(str_, pattern = " ")[[1]][6] # you will have to adjust the second number here if paths change
}

design = read.csv(file.path(IMPACD_R_folder, "annotations_cleaned.csv"))
design = merge(design, fcsFiles, by = "sampleName")

# Reduce to the samples that we have to batch 1
design = dplyr::filter(design, batch == 1, sample_id != "IC")

# ------------------------------------------------------------------------------
# Get gate definitions and process them
# ------------------------------------------------------------------------------
gate_definitions = read_excel("../KhorDS_CD4Tcells_Marker_Thresholds_Per_Batch_modified.xlsx", sheet = 'Final Thresholds')

batch_number = paste0("Batch", "_", "1")

gate_definitions = dplyr::select(gate_definitions, Annotation, Name, "Marker Type", `Gate Position` = all_of(batch_number))

# Save data
usethis::use_data(gate_definitions)

processed_gates = processGates(gate_definitions, transList_wsp)

gate_metadata = processed_gates[[1]]
gates = processed_gates[[2]]


# ------------------------------------------------------------------------------
# Transform and save fcs files
# ------------------------------------------------------------------------------
fcsDir <- system.file("extdata",
                       package="IMPACD")

design$sampleName = paste0("sample_", rownames(design))
design$fcsName = paste0("sample_", rownames(design), ".fcs")

# Writing the raw fcs file to the IMPACD extdata directory,
for (iRow in 1:nrow(design)) {
    ff = read.FCS(design$fileName[iRow])
    write.FCS(ff, file = file.path(fcsDir, design$fcsName[iRow]))
}


# Test that we can load and transform the fcs files
fcsfiles <- list.files(fcsDir, full = TRUE, pattern = "*.fcs")

fs = read.flowSet(fcsfiles)

for (fcs_ in names(fs@frames)) {
  print(fcs_)
  fs@frames[[fcs_]] = flowCore::transform(fs@frames[[fcs_]], transList_wsp)
}

# Verify everything looks good
gs = GatingSet(fs)

# ------------------------------------------------------------------------------
# Save annotation
# ------------------------------------------------------------------------------
annotation = dplyr::select(design, sampleName, fcsName, condition, age, gender)

usethis::use_data(annotation)

devtools::document()
