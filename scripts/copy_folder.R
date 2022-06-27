
# Copy folder -------------------------------------------------------------

# set directory
directory = "/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/"
setwd(directory)

# Transfer all files  -----------------------------------------------------

# out
system("gsutil cp rds_objects/gse* gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/rds_objects/")
# in (be careful - will overwrite current files)
system("gsutil cp -r gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/input_data/doublet_ID_6_22.tsv input_data/")


# Transfer specific folders  ----------------------------------------------

# helper functions
transfer_folder_out = function(folder) {
  system(paste0("gsutil cp -r ", folder, " gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/"))
}

transfer_folder_in = function(folder) {
  system(paste0("gsutil cp -r gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/", folder, " ."))
}

# execution
transfer_folder_out("scripts")
transfer_folder_out("figures")
transfer_folder_out("tables")

