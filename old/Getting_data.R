# data management libraries
library(readxl)
library(stringr)
library(dplyr)
# genetic libraries
library(bold)
library(seqinr)
library(rentrez)

species <- read_xlsx("data/species.xlsx")

# reading data from soortenlijst (https://doi.org/10.5281/zenodo.5596205)
# Duistermaat, H., Sparrius, L.B., & Denters, T. (2021). Standaardlijst van de Nederlandse Flora 2020 [Data set]. In Gorteria. Zenodo. https://doi.org/10.5281/zenodo.5596206
soortenlijst_nl <- read_xlsx("data/SL2020 Checklist Flora NL.xlsx", sheet = 2)
# splitting genus and species
soortenlijst_nl$genus <- word(soortenlijst_nl$`Wetenschappelijke naam`,1)
soortenlijst_nl$species <- str_replace(soortenlijst_nl$`Wetenschappelijke naam`, paste0(soortenlijst_nl$genus, " "), "")
# get data from bold on marker
bold_records <- bold_seq(marker = "trnL")
# join soortenlijst with bold-data
nl_soorten_bold <- soortenlijst_nl %>%
     select(soortnummer,`Wetenschappelijke naam`,`Nederlandse naam`, Familie, genus, species) %>%
     inner_join(bold_records, by = c("Wetenschappelijke naam" = "identification"))


#write.fasta(nl_soorten_bold,"nl_soorten_bold.fasta")



# get genbank data
entrez_db_summary("sra")

# check: https://fkeck.github.io/refdb/articles/ncbi_bold.html

