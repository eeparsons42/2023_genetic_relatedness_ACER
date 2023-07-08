install.packages('adegenet')
install.packages('ape')
install.packages('data.table')
install.packages('dbplyr')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('knitr')
install.packages('maps')
install.packages('mapproj')
install.packages('optparse')
install.packages('poppr')
install.packages('RColorBrewer')
install.packages('RPostgres')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

install.packages('tidyr')
install.packages('vcfR')
install.packages('vegan')
install.packages('yarrr')



suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dbplyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("maps"))
suppressPackageStartupMessages(library("mapproj"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("RPostgres"))
suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("yarrr"))
theme_set(theme_bw())


DEFAULT_MISSING_NUMERIC_VALUE <- -9.000000;

option_list <- list(
  make_option(c("--database_connection_string"), action="store", dest="database_connection_string", help="Corals (stag) database connection string"),
  make_option(c("--input_affy_metadata"), action="store", dest="input_affy_metadata", help="Affymetrix 96 well plate input file"),
  make_option(c("--input_pop_info"), action="store", dest="input_pop_info", help="Population information input file"),
  make_option(c("--input_vcf"), action="store", dest="input_vcf", help="VCF input file"),
  make_option(c("--output_nj_phylogeny_tree"), action="store", dest="output_nj_phylogeny_tree", default=NULL, help="Flag to plot neighbor-joining phylogeny tree"),
  make_option(c("--output_stag_db_report"), action="store", dest="output_stag_db_report", help="Flag to output stag db report file")
)

setwd("C:/Users/parsonsee/Desktop/Genetic analysis practice")
getwd()

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_file_path = function(dir, dup_removed_DB_snps.vcf.gz) {
  file_path = paste(dir, dup_removed_DB_snps.vcf.gz, sep="/");
  return(file_path);
}

get_database_connection <- function(db_conn_string) {
  # Instantiate database connection.
  # The connection string has this format:
  # postgresql://user:password@host/dbname
  conn_items <- strsplit(db_conn_string, "://")[[1]];
  string_needed <- conn_items[2];
  items_needed <- strsplit(string_needed, "@")[[1]];
  user_pass_string <- items_needed[1];
  host_dbname_string <- items_needed[2];
  user_pass_items <- strsplit(user_pass_string, ":")[[1]];
  host_dbname_items <- strsplit(host_dbname_string, "/")[[1]];
  user <- user_pass_items[1];
  pass <- user_pass_items[2];
  host <- host_dbname_items[1];
  dbname <- host_dbname_items[2];
  conn <- DBI::dbConnect(RPostgres::Postgres(), host=host, port="5432", dbname=dbname, user=user, password=pass);
  return (conn);
}

log_data_frame <- function(name, df) {
  cat("\n", name, ":\n");
  show(df);
  cat("\n\n");
}

time_elapsed <- function(start_time) {
  cat("Elapsed time: ", proc.time() - start_time, "\n\n");
}  

time_elapsed <- function(start_time) {
  cat("Elapsed time: ", proc.time() - start_time, "\n\n");
}

time_start <- function(msg) {
  start_time <- proc.time();
  cat(msg, "...\n");
  return(start_time);
}

write_data_frame <- function(dir, dup_removed_DB_snps.vcf.gz, data_frame) {
  cat("\nWriting file: ", dup_removed_DB_snps.vcf.gz, "\n");
  file_path <- get_file_path(dir, dup_removed_DB_snps.vcf.gz);
  write.table(data_frame, file=file_path, quote=FALSE, row.names=FALSE, sep="\t");
}


# Prepare for processing.
output_data_dir = "output_data_dir";
output_plots_dir = "output_plots_dir";
# Read in VCF input file.
start_time <- time_start("Reading VCF input");
vcf <- read.vcfR(opt$input_vcf);
time_elapsed(start_time);


# Convert VCF file into a genind for the Poppr package.
start_time <- time_start("Converting VCF data to a genind object");
genind_obj <- vcfR2genind(vcf);
cat("\ngenind_obj:\n");
genind_obj
cat("\n\n");
time_elapsed(start_time);


# Add population information to the genind object.
population_info_data_table <- read.table(opt$popInfo_p1, check.names=FALSE, header=F, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t", quote="");
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region");
cat("\npopulation_info_data_table:\n");
population_info_data_table
cat("\n\n");
#write_data_frame(output_data_dir, "population_info_data_table", population_info_data_table);
genind_obj@pop <- as.factor(population_info_data_table$region);
strata(genind_obj) <- data.frame(pop(genind_obj));

start_time <- time_start("Converting the genind object to a genclone object");
genind_clone <- as.genclone(genind_obj);
cat("\ngenind_clone:\n");
genind_clone
cat("\n\n");
time_elapsed(start_time);

# Calculate the bitwise distance between individuals.
start_time <- time_start("Calculating the bitwise distance between individuals");
bitwise_distance <- bitwise.dist(genind_clone);
time_elapsed(start_time);

# Multilocus genotypes (threshold of 3.2%).
cat("\nFiltering multilocus genotypes with threshold of 3.2%...\n\n");
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032;

# Create list of MLGs.
cat("\nCreating list of mlg_ids...\n\n");
mlg_ids <- mlg.id(genind_clone);

# Read user's Affymetrix 96 well plate tabular file.
affy_metadata_data_frame <- read.table(opt$Plate1_genotypingReport.Rmd, header=FALSE, stringsAsFactors=FALSE, sep="\t", na.strings=c("", "NA"), quote="");
colnames(affy_metadata_data_frame) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                                        "region", "latitude", "longitude", "geographic_origin", "colony_location",
                                        "depth", "disease_resist", "bleach_resist", "mortality","tle",
                                        "spawning", "collector_last_name", "collector_first_name", "organization", "collection_date",
                                        "email", "seq_facility", "array_version", "public", "public_after_date",
                                        "sperm_motility", "healing_time", "dna_extraction_method", "dna_concentration", "registry_id",
                                        "result_folder_name", "plate_barcode");
affy_metadata_data_frame$user_specimen_id <- as.character(affy_metadata_data_frame$user_specimen_id);
log_data_frame("affy_metadata_data_frame", affy_metadata_data_frame);
user_specimen_ids <- as.character(affy_metadata_data_frame$user_specimen_id);
cat("\nuser_specimen_ids:\n", toString(user_specimen_ids), "\n\n");
# The specimen_id_field_call_data_table looks like this:
# user_specimen_ids V2
# 1090              prolifera
# 1091              prolifera
cat("\nCreating specimen_id_field_call_data_table...\n");
specimen_id_field_call_data_table <- data.table(user_specimen_ids, affy_metadata_data_frame$field_call);
# Rename the user_specimen_ids column.
setnames(specimen_id_field_call_data_table, c("user_specimen_ids"), c("user_specimen_id"));
# Rename the V2 column.
setnames(specimen_id_field_call_data_table, c("V2"), c("field_call"));