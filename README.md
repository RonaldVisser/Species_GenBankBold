# Species_GenBankBold

Getting data from GenBank and Bold for further analyses

# Usage

The file [getting_combining_data.R](https://github.com/RonaldVisser/Species_GenBankBold/blob/main/getting_combining_data.R) is used to get data from GenBank and Bold on basis os the trnL marker. The species list of plants (see data) is used to select only plant DNA.

[ncbi_download.R](https://github.com/RonaldVisser/Species_GenBankBold/blob/main/ncbi_download.R) is a function to download data from the NCBI databank, used by the aforementioned script.

[refdb_import_BOLD.R](https://github.com/RonaldVisser/Species_GenBankBold/blob/main/refdb_import_BOLD.R) is an adapted version of the function from <https://github.com/fkeck/refdb/blob/main/R/import_NCBI.R>. See also:

Keck, F., & Altermatt, F. (2023). Management of DNA reference libraries for barcoding and metabarcoding studies with the R package refdb. *Molecular Ecology Resources*, 23(2), 511-518. <https://doi.org/10.1111/1755-0998.13723>.

# Data

The file SL2020 Checklist Flora NL.xlsx was taken from

Duistermaat, H., Sparrius, L.B., & Denters, T. (2021). *Standaardlijst van de Nederlandse Flora 2020* [Data set]. In Gorteria. Zenodo. <https://doi.org/10.5281/zenodo.5596206>
