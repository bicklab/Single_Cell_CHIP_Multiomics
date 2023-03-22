
# Copy folder -------------------------------------------------------------

# set directory
directory = "~/Alyssa/Single_Cell_CHIP/"
setwd(directory)

# Transfer all files  -----------------------------------------------------

# gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46 is bicklab-aps
# gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d is Alyssa_Parker
# out
system("gsutil cp rds_objects/gse* gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/rds_objects/")
# in (be careful - will overwrite current files)
system("gsutil cp -r gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/Single_Cell_CHIP/input_data/doublet_ID_6_22.tsv input_data/")

system("gsutil cp  gs://bicklab-main-storage/Workspaces/bicklab-chip-scRNAseq/integrated_seurat.rds ../Single_Cell_CHIP/rds_objects/")
system("gsutil cp -r ../* gs://fc-e44e4e80-489e-4c0f-a2a7-a96718d54f46/")
system("gsutil cp gs://bicklab-main-storage/Users/Brett_Heimlich/shared/seurat_DF_Harmony_KEEP.rds rds_objects/")
system("gsutil cp gs://bicklab-main-storage/Users/Brett_Heimlich/shared/seurat_DF.rds rds_objects/")

system("gsutil cp rds_objects/seurat_df_harmony_keep_0_75_annotated.rds gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/rds_objects/")


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

