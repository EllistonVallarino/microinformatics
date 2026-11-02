#' @title Microinformatics
#' @description
#' Package for cleaning and reshaping MOTHUR outputs for easier analysis and visualization.
#'
#' @param clean_shared, MOTHUR output .shared, without the OTU/ASVRare1, as an object
#' @param taxonomy_file, MOTHUR output taxonomy file .cons.taxonomy as an object
#'
#' @export
get_taxclean <- function(clean_shared, taxonomy_file){

  taxonomy_cons <- taxonomy_file
  clean_asvs <- colnames(clean_shared)
  clean_asvs <- clean_asvs[4:((length(clean_asvs)))]

  taxonomy_cleaned <- filter(taxonomy_cons, OTU %in% clean_asvs)

  write_tsv(taxonomy_cleaned, file = 'final_season_space.asv.ASV.cons.cleaned.taxonomy')

  taxonomy_cleaned <-  taxonomy_cleaned %>%
    rename_all(tolower)
  return(taxonomy_cleaned)

}
