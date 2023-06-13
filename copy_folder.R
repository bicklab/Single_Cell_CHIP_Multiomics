
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
system("gsutil cp input_data/pANN_frame.tsv gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/input_data/")

system("gsutil cp gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/rds_objects/DNMT3Amutant_seu.rds rds_objects/ ")
system("gsutil cp gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/rds_objects/TET2mutant_seu.rds rds_objects/ ")

system("gsutil cp gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/h5ad_files/*_for_seurat.h5ad h5ad_files/")

system("gsutil cp scripts/* gs://fc-secure-0a0db300-cd6b-49a0-a17e-633b998bfc65/CH_Multiomics/Code_folder/")
system("gsutil cp -r tables/meta_pseudobulk gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/tables/")
system("gsutil cp -r demultiplexed/* gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/demultiplexed/")

system("gsutil cp -r figures gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/")
system("gsutil cp gs://fc-secure-084dba6b-e7ee-4ebc-9fad-59c87bc8d1cc/Maester_data_PB/seurat_with_clones_5.11.23.rds rds_objects/")

system("gsutil ls gs://fc-secure-0a0db300-cd6b-49a0-a17e-633b998bfc65/CH_Multiomics/Code_folder/")
system("gsutil cp -r tables gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/")


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

