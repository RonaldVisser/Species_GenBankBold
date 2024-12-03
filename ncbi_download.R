#
#' Function for downloading ncbi-data
#' The function ncbi_download() was created to loop through a list of species and find data in the ncbi-database.
#' Dowloading data using refdb::refdb_import_NCBI() often results in errors in the connection. Therefore errors are handled by breaks.
#' The default break is an error occurs is 60 seconds, after that the function will try again to download data.
#' De default wait between two download requests is set to 10 seconds. This can be adapted if needed. The prefect value can be found using the famous trial and error method.
#' The function runs through each row of the dataset and uses the full botanical name of the plant to find the matching data string and with the search_string this can be expanded
#'
#' @param species_list vector with list of species, using the full botanical name (Latin)
#' @param search_string custom search string to combine with `species_list`. Default is null
#' @param error_wait numeric value setting the number of seconds to wait when downloading data from ncbi results in an error. Default value = 60
#' @param wait_time numeric value setting the number of seconds to wait between two download rounds. Default value  = 10
#' @returns ncbi_data tibble with ncbi_data for species list
#' @author Ronald M. Visser
#' @examples
#'
#'
#'
ncbi_download <- function(species_list, search_string = NULL, error_wait = 60, wait_time = 10) {
     for (r in 1:length(species_list)) {
          if (is.null(species_list) | !is.vector(species_list)) {
               print("You need to supply a species list as a vector")
          }
          message(paste0("Finding data for ", species_list[r]))
          if (!is.null(search_string)){
               search_data <- paste0(species_list[r], " ", search_string)
          } else {
               search_data <- species_list[r]
               }
          tryCatch(
               {ncbi_temp <- refdb_import_NCBI(search_data)},
               error = function(e) {
                    # add break of 60 seconds in case of error and retry
                    Sys.sleep(error_wait)
                    ncbi_temp <- refdb_import_NCBI(search_data)})
          if (!is.null(nrow(ncbi_temp))) {
               if (exists("ncbi_data")) {
                    ncbi_data <- rbind(ncbi_data, ncbi_temp)
               } else {
                    ncbi_data <- ncbi_temp
               }
          }
          Sys.sleep(wait_time) # to prevent errors: data usage is limited to max 3 queries per second: https://www.ncbi.nlm.nih.gov/home/about/policies/
     }
     if (exists("ncbi_data")) {
          ncbi_data }
     else {
          message("No data was found")
     }
}

# trnL_ncbi_temp <- refdb_import_NCBI(paste0("trnL voucher ", data_set$name[i]))
