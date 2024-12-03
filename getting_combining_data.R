### Downloading and combining data from NCBI Genbank and BOLD
### Author: Ronald Visser with input from Liam Oskam and Anja Fischer
###
# data management libraries
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
# gene data
library(refdb)
library(ape)
library(bold)
library(seqinr)
library(rentrez)
# special functions
source("ncbi_download.R") # function to donwnload from ncbi
source("refdb_import_BOLD.R") # adapted version of import db to include trnL gene

# reading data from soortenlijst (https://doi.org/10.5281/zenodo.5596205)
# Duistermaat, H., Sparrius, L.B., & Denters, T. (2021). Standaardlijst van de Nederlandse Flora 2020 [Data set]. In Gorteria. Zenodo. https://doi.org/10.5281/zenodo.5596206
soortenlijst_nl <- read_xlsx("data/SL2020 Checklist Flora NL.xlsx", sheet = 2)
# splitting genus and species
soortenlijst_nl$genus <- word(soortenlijst_nl$`Wetenschappelijke naam`,1)
soortenlijst_nl$species <- str_replace(soortenlijst_nl$`Wetenschappelijke naam`, paste0(soortenlijst_nl$genus, " "), "")
# limit columns
soortenlijst <- soortenlijst_nl %>% select(`Wetenschappelijke naam`, genus, species)
colnames(soortenlijst) <- c("name", "genus", "species")
soortenlijst$database <- 'floron'
# save to RDS for easy loading
saveRDS(soortenlijst, "data/soortenlijst.rds")
# load data
soortenlijst <- readRDS("data/soortenlijst.rds")


# data from RADAR:
# reference: Otto Brinkkemper: https://www.cultureelerfgoed.nl/onderwerpen/bronnen-en-kaarten/overzicht/zadendatabase
radar_species <- read_xlsx("data/species_in_RADAR_20211209.xlsx")
# removing problematic names
radar_species <- radar_species %>% filter(Name != 'Prunus domestica subsp. insititia (Gro-6|12)')
radar_species <- radar_species %>% filter(Name != 'Vaccinium (ex. Oxycoccus|Vaccinium)')
# string split on |
#radar_species$Name[str_detect(radar_species$Name, "\\|")]
radar_species <- unique(separate_longer_delim(radar_species, Name, delim= "|"))
radar_species$genus <- word(radar_species$Name,1)
radar_species$species <- str_replace(radar_species$Name, paste0(radar_species$genus, " "), "")
radar_species$species[radar_species$genus==radar_species$species] <- NA
radar_species <- radar_species %>% filter(!is.na(species))
colnames(radar_species) <- c("name", "genus", "species")
radar_species$database <- 'radar'
# save to RDS for easy loading
saveRDS(radar_species, "data/radar_species.rds")
# load data
radar_species <- readRDS("data/radar_species.rds")

species_combined <- unique(rbind(soortenlijst, radar_species))

species_combined <- species_combined %>%
     mutate(n = TRUE) %>%
     pivot_wider(names_from = database, values_from = n, values_fill = FALSE)

# save to RDS for easy loading
saveRDS(species_combined, "data/species_combined.rds")
write.csv(species_combined, "data/species_combinded.csv")
# load data (if needed)
species_combined <- readRDS("data/species_combined.rds")


# Getting data from NCBI --------------------------------------------------

# download from ncbi often gives HTTP failure: 400 and timeout errors, therefore the process has to be repeated several times to make sure that all data is donwloaded
# other errors: Could not resolve host: eutils.ncbi.nlm.nih.gov

split_parts <- 16 # number of parts to split the list of species in
species_combined_split <- split(species_combined, factor(sort(rank(row.names(species_combined))%% split_parts)))

# create dir for log files
dir.create(paste0("log/", format(Sys.Date(), "%Y%m%d")))

#' error handling is partly by waiting 60 seconds on error of fetching data from ncbi (added 20241130)
#' input is dataset with 5 variables: name (full botanical name), genus, species, floron (present in floron list), radar (present in radar)
#' the function runs through each row of the dataset and uses the full botanical name of the plant to find the matching trnL data string

# Round 7: unexplained 414 error with "Panicum barbipulvinatum", remove from split_part and downloaded separately,
#species_combined_split[[7]] <- subset(species_combined_split[[7]], name != "Panicum barbipulvinatum")
#trnL_ncbi_panicum <- refdb_import_NCBI("Panicum barbipulvinatum trnL voucher")
# Round 7: unexplained 414 error with "Ranunculus parviflorus", remove from split_part and downloaded separately,
#species_combined_split[[7]] <- subset(species_combined_split[[7]], name != "Ranunculus parviflorus")
#trnL_ncbi_ranunculus <- refdb_import_NCBI("Ranunculus parviflorus trnL voucher")

# Round 8: unexplained 414 error with "Prunus domestica subsp. insititia (Gro-11)"
#species_combined_split[[8]] <- subset(species_combined_split[[8]], name != "Prunus domestica subsp. insititia (Gro-11)")
#trnL_ncbi_prunus <- refdb_import_NCBI("Prunus domestica subsp. insititia (Gro-11) trnL voucher")

# loop through parts and download ncbi-data, process is logged in log folder with todays date
for (i in 1:split_parts) {
     log_file <- paste0("log/", format(Sys.Date(), "%Y%m%d"), "/ncbi_results_",  sprintf("%02d",i) ,".log")
     sink(log_file, append=TRUE, split=TRUE)
     #assign(paste0("trnL_ncbi", i), ncbi_download(species_combined_split[[toString(i-1)]]))
     assign(paste0("trnL_ncbi", sprintf("%02d",i)), ncbi_download(species_combined_split[[toString(i-1)]]$name, search_string = "trnL voucher", error_wait = 60, wait_time = 1))
     sink()
}
# if error, please call sink() and replace 1 in (i in 1:split_parts) to the number of the last run. Remove log file of last run


# to make sure all data is there, the data is combined and then filtered for unique records
ncbi_objects <- ls(pattern = "trnL_ncbi")
trnL_ncbi <- do.call("rbind", mget(ls(pattern = "trnL_ncbi")))
trnL_ncbi_filtered <- unique(trnL_ncbi)


# filter data on trnL gene and on species in list of species from radar and soortenlijst (extra check and leads to strong reduction)
trnL_ncbi_filtered <- trnL_ncbi_filtered %>%
     filter(gene=="trnL") %>%
     filter(species %in% species_combined$name)
# remove downloaded objects
rm(list=ncbi_objects)

# Reading data from BOLD --------------------------------------------------

# check presence of species in BOLD
species_in_bold <- bold_tax_name(name = species_combined$name)

species_in_bold_names <- species_in_bold$input[!is.na(species_in_bold$taxid)]

# read data from bold using refdb-package (slow)
for (i in 1:length(species_in_bold_names)) {
     message(paste0("Finding data for ", species_in_bold_names[i]))
     soorten_bold_temp <- refdb_import_BOLD(taxon=species_in_bold_names[i], marker="trnL")
     print(i)
     if (!is.null(nrow(soorten_bold_temp))) {
          if (exists("soorten_bold")) {
               soorten_bold <- rbind(soorten_bold, soorten_bold_temp)
          } else {
               soorten_bold <- soorten_bold_temp
          }
     }
}
# errors every now and then due to constraints in download

# filter marker and check with list of species from radar and soortenlijst (extra check)
trnL_bold <- soorten_bold %>%
     filter(markercode == "trnL") %>%
     filter(species %in% species_combined$name) %>%
     distinct()

# get data from bold on marker using bold-package, faster but less info
bold_records <- bold_seq(marker = "trnL")
# join soortenlijst with bold-data
nl_soorten_bold <- species_combined %>%
     select(name, genus, species) %>%
     inner_join(bold_records, by = c("name" = "identification"))

# Merging BOLD and NCBI ------------------------------------------------

trnL_bold_ncbi <- refdb_merge(trnL_bold, trnL_ncbi_filtered)

# Exporting data ----------------------------------------------------------


write.nexus.data(trnL_bold_ncbi, file="export/trnl_data.nexus",
                 format="protein")
write.dna(trnL_bold_ncbi, "export/trnl_data.dna")
write.csv(trnL_bold_ncbi, "export/trnl_data.csv", row.names = FALSE)

