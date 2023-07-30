### Downloading and combining data from NCBI Genbank and BOLD
### Author: Ronald Visser / Liam Oskam
###
# data management libraries
library(readxl)
library(stringr)
library(dplyr)
# gene data
library(refdb)
library(ape)
library(bold)
library(seqinr)
library(rentrez)

# reading data from soortenlijst (https://doi.org/10.5281/zenodo.5596205)
# Duistermaat, H., Sparrius, L.B., & Denters, T. (2021). Standaardlijst van de Nederlandse Flora 2020 [Data set]. In Gorteria. Zenodo. https://doi.org/10.5281/zenodo.5596206
soortenlijst_nl <- read_xlsx("data/SL2020 Checklist Flora NL.xlsx", sheet = 2)
# splitting genus and species
soortenlijst_nl$genus <- word(soortenlijst_nl$`Wetenschappelijke naam`,1)
soortenlijst_nl$species <- str_replace(soortenlijst_nl$`Wetenschappelijke naam`, paste0(soortenlijst_nl$genus, " "), "")

#trnL_ncbi <- refdb_import_NCBI("trnL voucher") # results in much data -> slow download

# Getting data from NCBI --------------------------------------------------

sink("log/ncbi_results.log", append=TRUE, split=TRUE)
for (i in 1:nrow(soortenlijst_nl)) {
     message(paste0("Finding data for ", soortenlijst_nl$`Wetenschappelijke naam`[i]))
     trnL_ncbi_temp <- refdb_import_NCBI(paste0("trnL voucher ", soortenlijst_nl$`Wetenschappelijke naam`[i]))
     if (!is.null(nrow(trnL_ncbi_temp))) {
          if (exists("trnL_ncbi")) {
               trnL_ncbi <- rbind(trnL_ncbi, trnL_ncbi_temp)
               } else {
                    trnL_ncbi <- trnL_ncbi_temp
               }
     }
     time.sleep(2) # to prevent errors (?)
}
sink()

trnL_ncbi1 <- trnL_ncbi # first round: 10499 observations (no log)
rm(trnL_ncbi)
# run lines 26-39 again
trnL_ncbi2 <- trnL_ncbi # second round 9045 observations
file.rename("log/ncbi_results.log", "log/ncbi_results_2.log")
rm(trnL_ncbi)
# run lines 26-39 again
trnL_ncbi3 <- trnL_ncbi # third round 6766 observations, with HTTP error
file.rename("log/ncbi_results.log", "log/ncbi_results_3.log")
rm(trnL_ncbi)
# run lines 26-39 again
trnL_ncbi4 <- trnL_ncbi # fourth round 10499 observations (same as first round)
file.rename("log/ncbi_results.log", "log/ncbi_results_4.log")
rm(trnL_ncbi)
# to make sure all data is there, the data is combined and then filtered for unique records
trnL_ncbi <- rbind(trnL_ncbi1, trnL_ncbi2, trnL_ncbi3, trnL_ncbi4)
trnL_ncbi <- unique(trnL_ncbi)
# filter data on trnL gene and on species in soortenlijst
trnL_ncbi <- trnL_ncbi %>%
     filter(gene=="trnL") %>%
     filter(species %in% soortenlijst_nl$`Wetenschappelijke naam`)



# Reading data from BOLD --------------------------------------------------

# read data from bold using refdb-package (slow)
for (i in 1:nrow(soortenlijst_nl)) {
     message(paste0("Finding data for ", soortenlijst_nl$`Wetenschappelijke naam`[i]))
     soorten_bold_temp <- refdb_import_BOLD(taxon=soortenlijst_nl$`Wetenschappelijke naam`[i])
     if (!is.null(nrow(soorten_bold_temp))) {
          if (exists("soorten_bold")) {
               soorten_bold <- rbind(soorten_bold, soorten_bold_temp)
          } else {
               soorten_bold <- soorten_bold_temp
          }
     }
}
# filter marker and check with soortenlijst
trnL_bold <- soorten_bold %>%
     filter(markercode == "trnL") %>%
     filter(species %in% soortenlijst_nl$`Wetenschappelijke naam`) %>%
     distinct()

# get data from bold on marker using bold-package, faster but less info
bold_records <- bold_seq(marker = "trnL")
# join soortenlijst with bold-data
nl_soorten_bold <- soortenlijst_nl %>%
     select(soortnummer,`Wetenschappelijke naam`,`Nederlandse naam`, Familie, genus, species) %>%
     inner_join(bold_records, by = c("Wetenschappelijke naam" = "identification"))



# Merging BOLD and NCBI ------------------------------------------------

trnL_bold_ncbi <- refdb_merge(trnL_bold, trnL_ncbi)



# Exporting data ----------------------------------------------------------


write.nexus.data(trnL_bold_ncbi, file="export/trnl_data.nexus",
                 format="protein")
write.dna(trnL_bold_ncbi, "export/trnl_data.dna")
write.csv(trnL_bold_ncbi, "export/trnl_data.csv", row.names = FALSE)

