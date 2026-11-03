#' @title Microinformatics
#Last Edit: 2026-10-28
#This script contains various functions for the analysis of 16S and 18S data

#CLEANING TAXONOMY####
#' @title Microinformatics
#' @description
#' If the shared file was filtered, this function removes those ASVs from the taxonomy table.
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

#CLEANING SHARED FILE####
#' @title Microinformatics
#' @description
#' This function removes OTU/ASVRare1 from the table
#'
#'
#' @param shared_file, MOTHUR output shared file
#' @export
get_shared <-  function(shared_file){

  clean_shared <- shared_file
  names(clean_shared)[length(clean_shared)] <- 'ASVRare1'

  return(clean_shared)
}

#get alpha####
#' @title Microinformatics
#' @description
#' This function selects only the diversity metrics from the MOTHUR file ending in .summary, product of the MOTHUR function summary.single.
#' It recovers: group, sobs, shannon, invsimpson, and coverage.
#' @param alpha_file, The MOTHUR .summary file, output of the summary.single function.
#' @export
get_alpha <- function(alpha_file){

  alpha <- alpha_file%>%
    select(., group, sobs, shannon, invsimpson, coverage)

  return(alpha)

}

#SEPARATING TAXONOMY FILE USING THE 6-LEVEL MODEL####
##genus####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For genus
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_genus_silva <- function(taxclean){
  taxa_filtered <- taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kindom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, genus) %>%
    mutate(genus = str_replace(genus, '\\(100\\)', ''))
  return(taxa_filtered)
}

##family####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For family
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_family_silva <- function(taxclean){
  taxa_filtered <- taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kindom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, family) %>%
    mutate(family = str_replace(family, '\\(100\\)', ''))
  return(taxa_filtered)
}

##class####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For class
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_class_silva <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kindom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, class) %>%
    mutate(class = str_replace(class, '\\(100\\)', ''))


  return(taxa_filtered)
}

##order####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For order
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_order_silva <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kindom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, order) %>%
    mutate(order = str_replace(order, '\\(100\\)', ''))


  return(taxa_filtered)
}



##phylum####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For phylum
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_phylum_silva <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kindom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, phylum) %>%
    mutate(phylum = str_replace(phylum, '\\(100\\)', ''))


  return(taxa_filtered)
}


##kingdom####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 6-level taxonomy model (i.e. SILVA).
#' For kingdom
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_kingdom <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('kingdom', 'phylum', 'class', 'order', 'family',
                                'genus'), sep = ';') %>%
    select(otu, size, kingdom) %>%
    mutate(kingdom = str_replace(kingdom, '\\(100\\)', ''))


  return(taxa_filtered)
}

#separating taxonomy file for Eukaryotes using the 9-level taxonomy model####
##species####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For species
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_species_pr2 <- function(taxclean){
  taxa_filtered <- taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, species) %>%
    mutate(genus = str_replace(species, '\\(100\\)', ''))
  return(taxa_filtered)
}

##genus####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For genus
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_genus_pr2 <- function(taxclean){
  taxa_filtered <- taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, genus) %>%
    mutate(genus = str_replace(genus, '\\(100\\)', ''))
  return(taxa_filtered)
}


##family####
#' @title Microinformatics
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For family
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_family_pr2 <- function(taxclean){
  taxa_filtered <- taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, family) %>%
    mutate(family = str_replace(family, '\\(100\\)', ''))
  return(taxa_filtered)
}

##order####
#' @title Microinformatics #' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For order
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_order <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, order) %>%
    mutate(order = str_replace(order, '\\(100\\)', ''))


  return(taxa_filtered)
}


##class####
#' @title Microinformatics #'
#' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For class
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_class_pr2 <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, class) %>%
    mutate(class = str_replace(class, '\\(100\\)', ''))


  return(taxa_filtered)
}
##subdivision####
#' @title Microinformatics #' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For division
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_subdivision_pr2 <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, subdivision) %>%
    mutate(subdivision = str_replace(subdivision, '\\(100\\)', ''))


  return(taxa_filtered)
}

##division####
#' @title Microinformatics #' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For division
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_supergroup_pr2 <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, division) %>%
    mutate(division = str_replace(division, '\\(100\\)', ''))


  return(taxa_filtered)
}

##supergroup####
#' @title Microinformatics #' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For supergroup
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_supergroup_pr2 <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, supergroup) %>%
    mutate(supergroup = str_replace(supergroup, '\\(100\\)', ''))


  return(taxa_filtered)
}

##domain####
#' @title Microinformatics #' @description
#' This function retains one specific taxonomic level from the taxonomy file.
#' To be used with databases that use the 9-level taxonomy model (i.e. PR2).
#' For domain
#' @param taxclean, cleaned taxonomy table
#' @export
filter_taxa_domain_pr2 <- function(taxclean){
  taxa_filtered <-  taxclean%>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace(taxonomy, ';$', '')) %>%
    #slice_head() %>% pull(taxonomy) %>%
    separate(taxonomy, into = c('domain', 'supergroup','division', 'subdivision', 'class',
                                'order', 'family', 'genus', 'species'), sep = ';') %>%
    select(otu, size, domain) %>%
    mutate(domain = str_replace(domain, '\\(100\\)', ''))


  return(taxa_filtered)
}

#separating taxonomy file keeping OTU for databases using 6-level taxonomy####
##phylum####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 6-level taxonomy databases (i.e. SILVA).
#' For phylum
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_phylum_silva <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('kingdom', 'phylum', 'class', 'order', 'family',
                     'genus'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           phylum = str_replace(string = phylum,
                                pattern='(.*)',
                                replacement = '*\\1*'),
           phylum = str_replace(string=phylum,
                                pattern = '\\*(.*)_unclassified\\*',
                                replacement='Unclassified *\\1*'),
           taxon = glue('{phylum} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##class####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 6-level taxonomy databases (i.e. SILVA).
#' For class
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_class_silva <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('kingdom', 'phylum', 'class', 'order', 'family',
                     'genus'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           class = str_replace(string = class,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           class = str_replace(string=class,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{class} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##order####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 6-level taxonomy databases (i.e. SILVA).
#' For order
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_order_silva <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('kingdom', 'phylum', 'class', 'order', 'family',
                     'genus'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           order = str_replace(string = order,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           order = str_replace(string=order,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{order} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##family####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 6-level taxonomy databases (i.e. SILVA).
#' For family
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_family_silva <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('kingdom', 'phylum', 'class', 'order', 'family',
                     'genus'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           family = str_replace(string = family,
                                pattern='(.*)',
                                replacement = '*\\1*'),
           family = str_replace(string=family,
                                pattern = '\\*(.*)_unclassified\\*',
                                replacement='Unclassified *\\1*'),
           taxon = glue('{family} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##genus####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 6-level taxonomy databases (i.e. SILVA).
#' For genus
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_silva <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('kingdom', 'phylum', 'class', 'order', 'family',
                     'genus'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           genus = str_replace(string = genus,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           genus = str_replace(string=genus,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{genus} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

#separating taxonomy keeping OTU for databases using 9-level taxonomy model####
##subdivision####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For subdivision
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_subdivision_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           subdivision = str_replace(string = subdivision,
                                     pattern='(.*)',
                                     replacement = '*\\1*'),
           subdivision = str_replace(string=subdivision,
                                     pattern = '\\*(.*)_unclassified\\*',
                                     replacement='Unclassified *\\1*'),
           taxon = glue('{subdivision} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##division####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For division
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           division = str_replace(string = division,
                                  pattern='(.*)',
                                  replacement = '*\\1*'),
           division = str_replace(string=division,
                                  pattern = '\\*(.*)_unclassified\\*',
                                  replacement='Unclassified *\\1*'),
           taxon = glue('{division} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##domain####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For domain
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           domain = str_replace(string = domain,
                                pattern='(.*)',
                                replacement = '*\\1*'),
           domain = str_replace(string=domain,
                                pattern = '\\*(.*)_unclassified\\*',
                                replacement='Unclassified *\\1*'),
           taxon = glue('{domain} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)
}


##supergroup####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For supergroup
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           supergroup = str_replace(string = supergroup,
                                    pattern='(.*)',
                                    replacement = '*\\1*'),
           supergroup = str_replace(string=supergroup,
                                    pattern = '\\*(.*)_unclassified\\*',
                                    replacement='Unclassified *\\1*'),
           taxon = glue('{supergroup} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##class####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For class
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_class_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           class = str_replace(string = class,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           class = str_replace(string=class,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{class} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##order####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For order
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_order_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           order = str_replace(string = order,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           order = str_replace(string=order,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{order} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##family####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For family
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_family_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           family = str_replace(string = family,
                                pattern='(.*)',
                                replacement = '*\\1*'),
           family = str_replace(string=family,
                                pattern = '\\*(.*)_unclassified\\*',
                                replacement='Unclassified *\\1*'),
           taxon = glue('{family} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##genus####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For genus
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           genus = str_replace(string = genus,
                               pattern='(.*)',
                               replacement = '*\\1*'),
           genus = str_replace(string=genus,
                               pattern = '\\*(.*)_unclassified\\*',
                               replacement='Unclassified *\\1*'),
           taxon = glue('{genus} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)

}

##species####
#' @title Microinformatics #' @description
#' This function separate the taxonomy file, keeping the OTU number along with the taxa.
#' This is the function version to be used with 9-level taxonomy databases (i.e. PR2).
#' For species
#' @param taxclean, cleaned taxonomy table
#' @export
taxa_otu_genus_pr2 <- function(taxclean){

  taxa_clean <- taxonomy_cleaned %>%
    select("otu", 'taxonomy') %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ';$','')) %>%
    separate(taxonomy,
             into= c('domain', 'supergroup','division', 'subdivision', 'class',
                     'order', 'family', 'genus', 'species'), sep = ';') %>%
    mutate(pretty_otu=str_replace(string=otu,
                                  pattern='tu0*',
                                  replacement = 'TU'),
           species = str_replace(string = species,
                                pattern='(.*)',
                                replacement = '*\\1*'),
           species = str_replace(string=species,
                                pattern = '\\*(.*)_unclassified\\*',
                                replacement='Unclassified *\\1*'),
           taxon = glue('{species} ({pretty_otu})')) %>%
    select(otu, taxon)

  return(taxa_clean)
}

#get transposed shared####
#' @title Microinformatics #' @description
#' This function transposes the shared file so that each column is a sample and each row is a an ASV
#'
#' @param shared, MOTHUR .shared table
#' @export
get_transposed_shared <- function(shared){
  t_shared <- transpose(shared)
  rownames(t_shared) <-  colnames(shared)
  colnames(t_shared) <-  rownames(shared)
  t_shared <- t_shared[c(-1,-3),]
  t_shared[,length(t_shared)+1] <- rownames(t_shared)
  colnames(t_shared) <- t_shared[1,]
  t_shared_dbl <- t_shared %>%
    sapply(., as.double) %>%
    data.frame(.)
  t_shared_dbl[,length(t_shared)] <- t_shared[,length(t_shared)]
  t_shared_dbl <- t_shared_dbl[-1,]
  return(t_shared_dbl)
}

#GET TOTAL COUNTS PER SAMPLE####
##subdivision####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For subdivision.
#' @param subdivision_per_sample, taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_subdivision <- function(subdivision_per_sample){
  subdivision_full_abundance <-  subdivision_per_sample[,4] %>%
    map_dbl(sum)

  subdivision_sample_percent <- data.frame(otu = subdivision_per_sample$otu,
                                           size = subdivision_per_sample$size,subdivision = subdivision_per_sample$subdivision)
  for(i in 4:length(subdivision_per_sample)){
    subdivision_sample_percent[i] <- (subdivision_per_sample[i])
  }


  subdivision_taxa_abundance <-  data.frame()
  unique_subdivision <- unique(subdivision_per_sample$subdivision)
  samples_subdivision <- colnames(subdivision_sample_percent[,4:length(subdivision_sample_percent)])
  for(i in 1:length(samples_subdivision)){
    for(j in 1:length(unique_subdivision)){
      subdivision_taxa_abundance[j,i] <- subdivision_sample_percent %>%
        filter(., subdivision == unique_subdivision[j]) %>%
        select(., samples_subdivision[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(subdivision_taxa_abundance) <- colnames(subdivision_sample_percent[,4:length(subdivision_sample_percent)])

  subdivision_abundance_per_sample <- data.frame(subdivision = unique_subdivision, subdivision_taxa_abundance)
  write_csv(subdivision_abundance_per_sample, file = './subdivision_abundance_per_sample.csv')
  subdivision_abundance_per_sample <- subdivision_abundance_per_sample %>%
    mutate(., subdivision = replace_na(subdivision, 'Rare'))

  return(subdivision_abundance_per_sample)

}

##class####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For class
#' @param class_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_class <- function(class_per_sample){
  class_full_abundance <-  class_per_sample[,4] %>%
    map_dbl(sum)

  class_sample_percent <- data.frame(otu = class_per_sample$otu,
                                     size = class_per_sample$size,class = class_per_sample$class)
  for(i in 4:length(class_per_sample)){
    class_sample_percent[i] <- (class_per_sample[i])
  }


  class_taxa_abundance <-  data.frame()
  unique_class <- unique(class_per_sample$class)
  samples_class <- colnames(class_sample_percent[,4:length(class_sample_percent)])
  for(i in 1:length(samples_class)){
    for(j in 1:length(unique_class)){
      class_taxa_abundance[j,i] <- class_sample_percent %>%
        filter(., class == unique_class[j]) %>%
        select(., samples_class[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(class_taxa_abundance) <- colnames(class_sample_percent[,4:length(class_sample_percent)])

  class_abundance_per_sample <- data.frame(class = unique_class, class_taxa_abundance)
  write_csv(class_abundance_per_sample, file = './class_abundance_per_sample.csv')
  class_abundance_per_sample <- class_abundance_per_sample %>%
    mutate(., class = replace_na(class, 'Rare'))

  return(class_abundance_per_sample)

}

##phylum####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For phylum
#' @param phylum_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_phylum <- function(phylum_per_sample){
  phylum_full_abundance <-  phylum_per_sample[,4] %>%
    map_dbl(sum)

  phylum_sample_percent <- data.frame(otu = phylum_per_sample$otu,
                                      size = phylum_per_sample$size,phylum = phylum_per_sample$phylum)
  for(i in 4:length(phylum_per_sample)){
    phylum_sample_percent[i] <- (phylum_per_sample[i])
  }


  phylum_taxa_abundance <-  data.frame()
  unique_phylum <- unique(phylum_per_sample$phylum)
  samples_phylum <- colnames(phylum_sample_percent[,4:length(phylum_sample_percent)])
  for(i in 1:length(samples_phylum)){
    for(j in 1:length(unique_phylum)){
      phylum_taxa_abundance[j,i] <- phylum_sample_percent %>%
        filter(., phylum == unique_phylum[j]) %>%
        select(., samples_phylum[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(phylum_taxa_abundance) <- colnames(phylum_sample_percent[,4:length(phylum_sample_percent)])

  phylum_abundance_per_sample <- data.frame(phylum = unique_phylum, phylum_taxa_abundance)
  write_csv(phylum_abundance_per_sample, file = './phylum_abundance_per_sample.csv')
  phylum_abundance_per_sample <- phylum_abundance_per_sample %>%
    mutate(., phylum = replace_na(phylum, 'Rare'))

  return(phylum_abundance_per_sample)

}

##division####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For division.
#' @param division_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_division <- function(division_per_sample){
  division_full_abundance <-  division_per_sample[,4] %>%
    map_dbl(sum)

  division_sample_percent <- data.frame(otu = division_per_sample$otu,
                                        size = division_per_sample$size,division = division_per_sample$division)
  for(i in 4:length(division_per_sample)){
    division_sample_percent[i] <- (division_per_sample[i])
  }


  division_taxa_abundance <-  data.frame()
  unique_division <- unique(division_per_sample$division)
  samples_division <- colnames(division_sample_percent[,4:length(division_sample_percent)])
  for(i in 1:length(samples_division)){
    for(j in 1:length(unique_division)){
      division_taxa_abundance[j,i] <- division_sample_percent %>%
        filter(., division == unique_division[j]) %>%
        select(., samples_division[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(division_taxa_abundance) <- colnames(division_sample_percent[,4:length(division_sample_percent)])

  division_abundance_per_sample <- data.frame(division = unique_division, division_taxa_abundance)
  write_csv(division_abundance_per_sample, file = './division_abundance_per_sample.csv')
  division_abundance_per_sample <- division_abundance_per_sample %>%
    mutate(., division = replace_na(division, 'Rare'))

  return(division_abundance_per_sample)

}


##supergroup####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For supergroup
#' @param supergroup_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_supergroup <- function(supergroup_per_sample){
  supergroup_full_abundance <-  supergroup_per_sample[,4] %>%
    map_dbl(sum)

  supergroup_sample_percent <- data.frame(otu = supergroup_per_sample$otu,
                                          size = supergroup_per_sample$size,supergroup = supergroup_per_sample$supergroup)
  for(i in 4:length(supergroup_per_sample)){
    supergroup_sample_percent[i] <- (supergroup_per_sample[i])
  }


  supergroup_taxa_abundance <-  data.frame()
  unique_supergroup <- unique(supergroup_per_sample$supergroup)
  samples_supergroup <- colnames(supergroup_sample_percent[,4:length(supergroup_sample_percent)])
  for(i in 1:length(samples_supergroup)){
    for(j in 1:length(unique_supergroup)){
      supergroup_taxa_abundance[j,i] <- supergroup_sample_percent %>%
        filter(., supergroup == unique_supergroup[j]) %>%
        select(., samples_supergroup[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(supergroup_taxa_abundance) <- colnames(supergroup_sample_percent[,4:length(supergroup_sample_percent)])

  supergroup_abundance_per_sample <- data.frame(supergroup = unique_supergroup, supergroup_taxa_abundance)
  write_csv(supergroup_abundance_per_sample, file = './supergroup_abundance_per_sample.csv')
  supergroup_abundance_per_sample <- supergroup_abundance_per_sample %>%
    mutate(., supergroup = replace_na(supergroup, 'Rare'))

  return(supergroup_abundance_per_sample)

}

##domain####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For domain
#' @param domain_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_domain <- function(domain_per_sample){
  domain_full_abundance <-  domain_per_sample[,4] %>%
    map_dbl(sum)

  domain_sample_percent <- data.frame(otu = domain_per_sample$otu,
                                      size = domain_per_sample$size,domain = domain_per_sample$domain)
  for(i in 4:length(domain_per_sample)){
    domain_sample_percent[i] <- (domain_per_sample[i])
  }


  domain_taxa_abundance <-  data.frame()
  unique_domain <- unique(domain_per_sample$domain)
  samples_domain <- colnames(domain_sample_percent[,4:length(domain_sample_percent)])
  for(i in 1:length(samples_domain)){
    for(j in 1:length(unique_domain)){
      domain_taxa_abundance[j,i] <- domain_sample_percent %>%
        filter(., domain == unique_domain[j]) %>%
        select(., samples_domain[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(domain_taxa_abundance) <- colnames(domain_sample_percent[,4:length(domain_sample_percent)])

  domain_abundance_per_sample <- data.frame(domain = unique_domain, domain_taxa_abundance)
  write_csv(domain_abundance_per_sample, file = './domain_abundance_per_sample.csv')
  domain_abundance_per_sample <- domain_abundance_per_sample %>%
    mutate(., domain = replace_na(domain, 'Rare'))

  return(domain_abundance_per_sample)

}


##family####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For family
#' @param family_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_family <- function(family_per_sample){
  family_full_abundance <-  family_per_sample[,4] %>%
    map_dbl(sum)

  family_sample_percent <- data.frame(otu = family_per_sample$otu,
                                      size = family_per_sample$size,family = family_per_sample$family)
  for(i in 4:length(family_per_sample)){
    family_sample_percent[i] <- (family_per_sample[i])
  }


  family_taxa_abundance <-  data.frame()
  unique_family <- unique(family_per_sample$family)
  samples_family <- colnames(family_sample_percent[,4:length(family_sample_percent)])
  for(i in 1:length(samples_family)){
    for(j in 1:length(unique_family)){
      family_taxa_abundance[j,i] <- family_sample_percent %>%
        filter(., family == unique_family[j]) %>%
        select(., samples_family[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(family_taxa_abundance) <- colnames(family_sample_percent[,4:length(family_sample_percent)])

  family_abundance_per_sample <- data.frame(family = unique_family, family_taxa_abundance)
  write_csv(family_abundance_per_sample, file = './family_abundance_per_sample.csv')
  family_abundance_per_sample <- family_abundance_per_sample %>%
    mutate(., family = replace_na(family, 'Rare'))

  return(family_abundance_per_sample)

}

##order####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For order
#' @param order_per_sample, taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_order <- function(order_per_sample){
  order_full_abundance <-  order_per_sample[,4] %>%
    map_dbl(sum)

  order_sample_percent <- data.frame(otu = order_per_sample$otu,
                                     size = order_per_sample$size,order = order_per_sample$order)
  for(i in 4:length(order_per_sample)){
    order_sample_percent[i] <- (order_per_sample[i])
  }


  order_taxa_abundance <-  data.frame()
  unique_order <- unique(order_per_sample$order)
  samples_order <- colnames(order_sample_percent[,4:length(order_sample_percent)])
  for(i in 1:length(samples_order)){
    for(j in 1:length(unique_order)){
      order_taxa_abundance[j,i] <- order_sample_percent %>%
        filter(., order == unique_order[j]) %>%
        select(., samples_order[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(order_taxa_abundance) <- colnames(order_sample_percent[,4:length(order_sample_percent)])

  order_abundance_per_sample <- data.frame(order = unique_order, order_taxa_abundance)
  write_csv(order_abundance_per_sample, file = './order_abundance_per_sample.csv')
  order_abundance_per_sample <- order_abundance_per_sample %>%
    mutate(., order = replace_na(order, 'Rare'))

  return(order_abundance_per_sample)

}



##genus####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For genus
#' @param genus_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_genus <- function(genus_per_sample){
  genus_full_abundance <-  genus_per_sample[,4] %>%
    map_dbl(sum)

  genus_sample_percent <- data.frame(otu = genus_per_sample$otu,
                                     size = genus_per_sample$size,genus = genus_per_sample$genus)
  for(i in 4:length(genus_per_sample)){
    genus_sample_percent[i] <- (genus_per_sample[i])
  }


  genus_taxa_abundance <-  data.frame()
  unique_genus <- unique(genus_per_sample$genus)
  samples_genus <- colnames(genus_sample_percent[,4:length(genus_sample_percent)])
  for(i in 1:length(samples_genus)){
    for(j in 1:length(unique_genus)){
      genus_taxa_abundance[j,i] <- genus_sample_percent %>%
        filter(., genus == unique_genus[j]) %>%
        select(., samples_genus[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(genus_taxa_abundance) <- colnames(genus_sample_percent[,4:length(genus_sample_percent)])

  genus_abundance_per_sample <- data.frame(genus = unique_genus, genus_taxa_abundance)
  write_csv(genus_abundance_per_sample, file = './genus_abundance_per_sample.csv')
  genus_abundance_per_sample <- genus_abundance_per_sample %>%
    mutate(., genus = replace_na(genus, 'Rare'))

  return(genus_abundance_per_sample)

}


##species####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For species
#' @param species_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_species <- function(species_per_sample){
  species_full_abundance <-  species_per_sample[,4] %>%
    map_dbl(sum)

  species_sample_percent <- data.frame(otu = species_per_sample$otu,
                                     size = species_per_sample$size,species = species_per_sample$species)
  for(i in 4:length(species_per_sample)){
    species_sample_percent[i] <- (species_per_sample[i])
  }


  species_taxa_abundance <-  data.frame()
  unique_species <- unique(species_per_sample$species)
  samples_species <- colnames(species_sample_percent[,4:length(species_sample_percent)])
  for(i in 1:length(samples_species)){
    for(j in 1:length(unique_species)){
      species_taxa_abundance[j,i] <- species_sample_percent %>%
        filter(., species == unique_species[j]) %>%
        select(., samples_species[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(species_taxa_abundance) <- colnames(species_sample_percent[,4:length(species_sample_percent)])

  species_abundance_per_sample <- data.frame(species = unique_species, species_taxa_abundance)
  write_csv(species_abundance_per_sample, file = './species_abundance_per_sample.csv')
  species_abundance_per_sample <- species_abundance_per_sample %>%
    mutate(., species = replace_na(species, 'Rare'))

  return(species_abundance_per_sample)

}

##ASVs 18s####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For ASVs 18S using phylum as the taxonomic level.
#' @param asv_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_per_sample_asv_18s <- function(asv_per_sample){
  asv_full_abundance <-  asv_per_sample[,4] %>%
    map_dbl(sum)

  asv_sample_percent <- data.frame(asv = asv_per_sample$otu,
                                   size = asv_per_sample$size,
                                   genus = asv_per_sample$phylum)
  for(i in 4:length(asv_per_sample)){
    asv_sample_percent[i] <- (asv_per_sample[i])
  }


  asv_taxa_abundance <-  data.frame()
  unique_asv <- unique(asv_per_sample$otu)
  samples_asv <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])
  for(i in 1:length(samples_asv)){
    for(j in 1:length(unique_asv)){
      asv_taxa_abundance[j,i] <- asv_sample_percent %>%
        filter(., asv == unique_asv[j]) %>%
        select(., samples_asv[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(asv_taxa_abundance) <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])

  asv_abundance_per_sample <- data.frame(asv = unique_asv, phylum=asv_per_sample$phylum,
                                         asv_taxa_abundance)
  write_csv(asv_abundance_per_sample, file = './asv_abundance_per_sample.csv')
  asv_abundance_per_sample <- asv_abundance_per_sample %>%
    mutate(., asv = replace_na(asv, 'Rare'))

  return(asv_abundance_per_sample)

}

##ASVs 16s####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For ASVs 16S using genus as the taxonomic level.
#' @param asv_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_per_sample_asv_16s <- function(asv_per_sample){
  asv_full_abundance <-  asv_per_sample[,4] %>%
    map_dbl(sum)

  asv_sample_percent <- data.frame(asv = asv_per_sample$otu,
                                   size = asv_per_sample$size,
                                   genus = asv_per_sample$genus)
  for(i in 4:length(asv_per_sample)){
    asv_sample_percent[i] <- (asv_per_sample[i])
  }


  asv_taxa_abundance <-  data.frame()
  unique_asv <- unique(asv_per_sample$otu)
  samples_asv <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])
  for(i in 1:length(samples_asv)){
    for(j in 1:length(unique_asv)){
      asv_taxa_abundance[j,i] <- asv_sample_percent %>%
        filter(., asv == unique_asv[j]) %>%
        select(., samples_asv[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(asv_taxa_abundance) <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])

  asv_abundance_per_sample <- data.frame(asv = unique_asv, genus=asv_per_sample$genus,
                                         asv_taxa_abundance)
  write_csv(asv_abundance_per_sample, file = './asv_abundance_per_sample.csv')
  asv_abundance_per_sample <- asv_abundance_per_sample %>%
    mutate(., asv = replace_na(asv, 'Rare'))

  return(asv_abundance_per_sample)

}

##For general taxon####
#' @title Microinformatics #' @description
#' The following functions are to calculate the total counts for each taxa in each sample.
#' The function is done for each taxonomic level.
#' For general use, no need to specify the taxonomic level.
#' @param taxon_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_total_count_per_sample_taxon <- function(taxon_per_sample){
  taxon_full_abundance <-  taxon_per_sample[,4] %>%
    map_dbl(sum)

  taxon_sample_percent <- data.frame(otu = taxon_per_sample$otu,
                                     size = taxon_per_sample$size,taxon = pull(taxon_per_sample[,3]))
  for(i in 4:length(taxon_per_sample)){
    taxon_sample_percent[i] <- (taxon_per_sample[i])
  }


  taxon_taxa_abundance <-  data.frame()
  unique_taxon <- unique(pull(taxon_per_sample[,3]))
  samples_taxon <- colnames(taxon_sample_percent[,4:length(taxon_sample_percent)])
  for(i in 1:length(samples_taxon)){
    for(j in 1:length(unique_taxon)){
      taxon_taxa_abundance[j,i] <- taxon_sample_percent %>%
        filter(., taxon == unique_taxon[j]) %>%
        select(., samples_taxon[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(taxon_taxa_abundance) <- colnames(taxon_sample_percent[,4:length(taxon_sample_percent)])

  taxon_abundance_per_sample <- data.frame(taxon = unique_taxon, taxon_taxa_abundance)
  write_csv(taxon_abundance_per_sample, file = './taxon_abundance_per_sample.csv')
  taxon_abundance_per_sample <- taxon_abundance_per_sample %>%
    mutate(., taxon = replace_na(taxon, 'Rare'))

  return(taxon_abundance_per_sample)

}

#GET RELATIVE ABUNDANCE PER SAMPLE####
##for general taxon####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For general use, no need to specify the taxonomic level.
#' @param taxon_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_taxon <- function(taxon_per_sample){
  taxon_full_abundance <-  taxon_per_sample[,4] %>%
    map_dbl(sum)

  taxon_sample_percent <- data.frame(otu = taxon_per_sample$otu,
                                          size = taxon_per_sample$size, taxon = pull(taxon_per_sample[,3]))
  for(i in 4:length(taxon_per_sample)){
    taxon_sample_percent[i] <- (taxon_per_sample[i]/taxon_full_abundance)*100
  }


  taxon_taxa_abundance <-  data.frame()
  unique_taxon <- unique(pull(taxon_per_sample[,3]))
  samples_taxon <- colnames(taxon_sample_percent[,4:length(taxon_sample_percent)])
  for(i in 1:length(samples_taxon)){
    for(j in 1:length(unique_taxon)){
      taxon_taxa_abundance[j,i] <- taxon_sample_percent %>%
        filter(., taxon == unique_taxon[j]) %>%
        select(., samples_taxon[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(taxon_taxa_abundance) <- colnames(taxon_sample_percent[,4:length(taxon_sample_percent)])

  taxon_abundance_per_sample <- data.frame(taxon = unique_taxon, taxon_taxa_abundance)
  write_csv(taxon_abundance_per_sample, file = './taxon_abundance_per_sample.csv')
  taxon_abundance_per_sample <- taxon_abundance_per_sample %>%
    mutate(., taxon = replace_na(taxon, 'Rare'))
  return(taxon_abundance_per_sample)

}

##supergroup####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For supergroup.
#' @param supergroup_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_supergroup <- function(supergroup_per_sample){
  supergroup_full_abundance <-  supergroup_per_sample[,4] %>%
    map_dbl(sum)

  supergroup_sample_percent <- data.frame(otu = supergroup_per_sample$otu,
                                          size = supergroup_per_sample$size,supergroup = supergroup_per_sample$supergroup)
  for(i in 4:length(supergroup_per_sample)){
    supergroup_sample_percent[i] <- (supergroup_per_sample[i]/supergroup_full_abundance)*100
  }


  supergroup_taxa_abundance <-  data.frame()
  unique_supergroup <- unique(supergroup_per_sample$supergroup)
  samples_supergroup <- colnames(supergroup_sample_percent[,4:length(supergroup_sample_percent)])
  for(i in 1:length(samples_supergroup)){
    for(j in 1:length(unique_supergroup)){
      supergroup_taxa_abundance[j,i] <- supergroup_sample_percent %>%
        filter(., supergroup == unique_supergroup[j]) %>%
        select(., samples_supergroup[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(supergroup_taxa_abundance) <- colnames(supergroup_sample_percent[,4:length(supergroup_sample_percent)])

  supergroup_abundance_per_sample <- data.frame(supergroup = unique_supergroup, supergroup_taxa_abundance)
  write_csv(supergroup_abundance_per_sample, file = './supergroup_abundance_per_sample.csv')
  supergroup_abundance_per_sample <- supergroup_abundance_per_sample %>%
    mutate(., supergroup = replace_na(supergroup, 'Rare'))
  return(supergroup_abundance_per_sample)

}

##domain####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For domain.
#' @param domain_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_domain <- function(domain_per_sample){
  domain_full_abundance <-  domain_per_sample[,4] %>%
    map_dbl(sum)

  domain_sample_percent <- data.frame(otu = domain_per_sample$otu,
                                      size = domain_per_sample$size,domain = domain_per_sample$domain)
  for(i in 4:length(domain_per_sample)){
    domain_sample_percent[i] <- (domain_per_sample[i]/domain_full_abundance)*100
  }


  domain_taxa_abundance <-  data.frame()
  unique_domain <- unique(domain_per_sample$domain)
  samples_domain <- colnames(domain_sample_percent[,4:length(domain_sample_percent)])
  for(i in 1:length(samples_domain)){
    for(j in 1:length(unique_domain)){
      domain_taxa_abundance[j,i] <- domain_sample_percent %>%
        filter(., domain == unique_domain[j]) %>%
        select(., samples_domain[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(domain_taxa_abundance) <- colnames(domain_sample_percent[,4:length(domain_sample_percent)])

  domain_abundance_per_sample <- data.frame(domain = unique_domain, domain_taxa_abundance)
  write_csv(domain_abundance_per_sample, file = './domain_abundance_per_sample.csv')
  domain_abundance_per_sample <- domain_abundance_per_sample %>%
    mutate(., domain = replace_na(domain, 'Rare'))
  return(domain_abundance_per_sample)

}

##division####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For division.
#' @param division_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export

get_rel_abun_per_sample_division <- function(division_per_sample){
  division_full_abundance <-  division_per_sample[,4] %>%
    map_dbl(sum)

  division_sample_percent <- data.frame(otu = division_per_sample$otu,
                                        size = division_per_sample$size,division = division_per_sample$division)
  for(i in 4:length(division_per_sample)){
    division_sample_percent[i] <- (division_per_sample[i]/division_full_abundance)*100
  }


  division_taxa_abundance <-  data.frame()
  unique_division <- unique(division_per_sample$division)
  samples_division <- colnames(division_sample_percent[,4:length(division_sample_percent)])
  for(i in 1:length(samples_division)){
    for(j in 1:length(unique_division)){
      division_taxa_abundance[j,i] <- division_sample_percent %>%
        filter(., division == unique_division[j]) %>%
        select(., samples_division[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(division_taxa_abundance) <- colnames(division_sample_percent[,4:length(division_sample_percent)])

  division_abundance_per_sample <- data.frame(division = unique_division, division_taxa_abundance)
  write_csv(division_abundance_per_sample, file = './division_abundance_per_sample.csv')
  division_abundance_per_sample <- division_abundance_per_sample %>%
    mutate(., division = replace_na(division, 'Rare'))
  return(division_abundance_per_sample)

}

##subdivision####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For subdivision.
#' @param subdivision_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_subdivision <- function(subdivision_per_sample){
  subdivision_full_abundance <-  subdivision_per_sample[,4] %>%
    map_dbl(sum)

  subdivision_sample_percent <- data.frame(otu = subdivision_per_sample$otu,
                                           size = subdivision_per_sample$size,subdivision = subdivision_per_sample$subdivision)
  for(i in 4:length(subdivision_per_sample)){
    subdivision_sample_percent[i] <- (subdivision_per_sample[i]/subdivision_full_abundance)*100
  }


  subdivision_taxa_abundance <-  data.frame()
  unique_subdivision <- unique(subdivision_per_sample$subdivision)
  samples_subdivision <- colnames(subdivision_sample_percent[,4:length(subdivision_sample_percent)])
  for(i in 1:length(samples_subdivision)){
    for(j in 1:length(unique_subdivision)){
      subdivision_taxa_abundance[j,i] <- subdivision_sample_percent %>%
        filter(., subdivision == unique_subdivision[j]) %>%
        select(., samples_subdivision[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(subdivision_taxa_abundance) <- colnames(subdivision_sample_percent[,4:length(subdivision_sample_percent)])

  subdivision_abundance_per_sample <- data.frame(subdivision = unique_subdivision, subdivision_taxa_abundance)
  write_csv(subdivision_abundance_per_sample, file = './subdivision_abundance_per_sample.csv')
  subdivision_abundance_per_sample <- subdivision_abundance_per_sample %>%
    mutate(., subdivision = replace_na(subdivision, 'Rare'))
  return(subdivision_abundance_per_sample)

}

##class####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For class.
#' @param class_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_class <- function(class_per_sample){
  class_full_abundance <-  class_per_sample[,4] %>%
    map_dbl(sum)

  class_sample_percent <- data.frame(otu = class_per_sample$otu,
                                     size = class_per_sample$size,class = class_per_sample$class)
  for(i in 4:length(class_per_sample)){
    class_sample_percent[i] <- (class_per_sample[i]/class_full_abundance)*100
  }


  class_taxa_abundance <-  data.frame()
  unique_class <- unique(class_per_sample$class)
  samples_class <- colnames(class_sample_percent[,4:length(class_sample_percent)])
  for(i in 1:length(samples_class)){
    for(j in 1:length(unique_class)){
      class_taxa_abundance[j,i] <- class_sample_percent %>%
        filter(., class == unique_class[j]) %>%
        select(., samples_class[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(class_taxa_abundance) <- colnames(class_sample_percent[,4:length(class_sample_percent)])

  class_abundance_per_sample <- data.frame(class = unique_class, class_taxa_abundance)
  write_csv(class_abundance_per_sample, file = './class_abundance_per_sample.csv')
  class_abundance_per_sample <- class_abundance_per_sample %>%
    mutate(., class = replace_na(class, 'Rare'))
  return(class_abundance_per_sample)

}


##order####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For order.
#' @param order_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_order <- function(order_per_sample){
  order_full_abundance <-  order_per_sample[,4] %>%
    map_dbl(sum)

  order_sample_percent <- data.frame(otu = order_per_sample$otu,
                                     size = order_per_sample$size,order = order_per_sample$order)
  for(i in 4:length(order_per_sample)){
    order_sample_percent[i] <- (order_per_sample[i]/order_full_abundance)*100
  }


  order_taxa_abundance <-  data.frame()
  unique_order <- unique(order_per_sample$order)
  samples_order <- colnames(order_sample_percent[,4:length(order_sample_percent)])
  for(i in 1:length(samples_order)){
    for(j in 1:length(unique_order)){
      order_taxa_abundance[j,i] <- order_sample_percent %>%
        filter(., order == unique_order[j]) %>%
        select(., samples_order[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(order_taxa_abundance) <- colnames(order_sample_percent[,4:length(order_sample_percent)])

  order_abundance_per_sample <- data.frame(order = unique_order, order_taxa_abundance)
  write_csv(order_abundance_per_sample, file = './order_abundance_per_sample.csv')
  order_abundance_per_sample <- order_abundance_per_sample %>%
    mutate(., order = replace_na(order, 'Rare'))

  return(order_abundance_per_sample)

}

##family####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For family.
#' @param family_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_family <- function(family_per_sample){
  family_full_abundance <-  family_per_sample[,4] %>%
    map_dbl(sum)

  family_sample_percent <- data.frame(otu = family_per_sample$otu,
                                      size = family_per_sample$size,family = family_per_sample$family)
  for(i in 4:length(family_per_sample)){
    family_sample_percent[i] <- (family_per_sample[i]/family_full_abundance)*100
  }


  family_taxa_abundance <-  data.frame()
  unique_family <- unique(family_per_sample$family)
  samples_family <- colnames(family_sample_percent[,4:length(family_sample_percent)])
  for(i in 1:length(samples_family)){
    for(j in 1:length(unique_family)){
      family_taxa_abundance[j,i] <- family_sample_percent %>%
        filter(., family == unique_family[j]) %>%
        select(., samples_family[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(family_taxa_abundance) <- colnames(family_sample_percent[,4:length(family_sample_percent)])

  family_abundance_per_sample <- data.frame(family = unique_family, family_taxa_abundance)
  write_csv(family_abundance_per_sample, file = './family_abundance_per_sample.csv')
  family_abundance_per_sample <- family_abundance_per_sample %>%
    mutate(., family = replace_na(family, 'Rare'))

  return(family_abundance_per_sample)

}


##genus####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For genus.
#' @param genus_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export

get_rel_abun_per_sample_genus <- function(genus_per_sample){
  genus_full_abundance <-  genus_per_sample[,4] %>%
    map_dbl(sum)

  genus_sample_percent <- data.frame(otu = genus_per_sample$otu,
                                     size = genus_per_sample$size,genus = genus_per_sample$genus)
  for(i in 4:length(genus_per_sample)){
    genus_sample_percent[i] <- (genus_per_sample[i]/genus_full_abundance)*100
  }


  genus_taxa_abundance <-  data.frame()
  unique_genus <- unique(genus_per_sample$genus)
  samples_genus <- colnames(genus_sample_percent[,4:length(genus_sample_percent)])
  for(i in 1:length(samples_genus)){
    for(j in 1:length(unique_genus)){
      genus_taxa_abundance[j,i] <- genus_sample_percent %>%
        filter(., genus == unique_genus[j]) %>%
        select(., samples_genus[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(genus_taxa_abundance) <- colnames(genus_sample_percent[,4:length(genus_sample_percent)])

  genus_abundance_per_sample <- data.frame(genus = unique_genus, genus_taxa_abundance)
  write_csv(genus_abundance_per_sample, file = './genus_abundance_per_sample.csv')
  genus_abundance_per_sample <- genus_abundance_per_sample %>%
    mutate(., genus = replace_na(genus, 'Rare'))

  return(genus_abundance_per_sample)

}

##species####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For species.
#' @param species_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export

get_rel_abun_per_sample_species <- function(species_per_sample){
  species_full_abundance <-  species_per_sample[,4] %>%
    map_dbl(sum)

  species_sample_percent <- data.frame(otu = species_per_sample$otu,
                                     size = species_per_sample$size,species = species_per_sample$species)
  for(i in 4:length(species_per_sample)){
    species_sample_percent[i] <- (species_per_sample[i]/species_full_abundance)*100
  }


  species_taxa_abundance <-  data.frame()
  unique_species <- unique(species_per_sample$species)
  samples_species <- colnames(species_sample_percent[,4:length(species_sample_percent)])
  for(i in 1:length(samples_species)){
    for(j in 1:length(unique_species)){
      species_taxa_abundance[j,i] <- species_sample_percent %>%
        filter(., species == unique_species[j]) %>%
        select(., samples_species[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(species_taxa_abundance) <- colnames(species_sample_percent[,4:length(species_sample_percent)])

  species_abundance_per_sample <- data.frame(species = unique_species, species_taxa_abundance)
  write_csv(species_abundance_per_sample, file = './species_abundance_per_sample.csv')
  species_abundance_per_sample <- species_abundance_per_sample %>%
    mutate(., species = replace_na(species, 'Rare'))

  return(species_abundance_per_sample)

}


##phylum####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For phylum
#' @param phylum_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_phylum <- function(phylum_per_sample){
  phylum_full_abundance <-  phylum_per_sample[,4] %>%
    map_dbl(sum)

  phylum_sample_percent <- data.frame(otu = phylum_per_sample$otu,
                                      size = phylum_per_sample$size,phylum = phylum_per_sample$phylum)
  for(i in 4:length(phylum_per_sample)){
    phylum_sample_percent[i] <- (phylum_per_sample[i]/phylum_full_abundance)*100
  }


  phylum_taxa_abundance <-  data.frame()
  unique_phylum <- unique(phylum_per_sample$phylum)
  samples_phylum <- colnames(phylum_sample_percent[,4:length(phylum_sample_percent)])
  for(i in 1:length(samples_phylum)){
    for(j in 1:length(unique_phylum)){
      phylum_taxa_abundance[j,i] <- phylum_sample_percent %>%
        filter(., phylum == unique_phylum[j]) %>%
        select(., samples_phylum[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(phylum_taxa_abundance) <- colnames(phylum_sample_percent[,4:length(phylum_sample_percent)])

  phylum_abundance_per_sample <- data.frame(phylum = unique_phylum, phylum_taxa_abundance)
  write_csv(phylum_abundance_per_sample, file = './phylum_abundance_per_sample.csv')
  phylum_abundance_per_sample <- phylum_abundance_per_sample %>%
    mutate(., phylum = replace_na(phylum, 'Rare'))

  return(phylum_abundance_per_sample)

}


##kingdom####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For kingdom.
#' @param kingdom_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_kingdom <- function(kingdom_per_sample){
  kingdom_full_abundance <-  kingdom_per_sample[,4] %>%
    map_dbl(sum)

  kingdom_sample_percent <- data.frame(otu = kingdom_per_sample$otu,
                                       size = kingdom_per_sample$size,kingdom = kingdom_per_sample$kingdom)
  for(i in 4:length(kingdom_per_sample)){
    kingdom_sample_percent[i] <- (kingdom_per_sample[i]/kingdom_full_abundance)*100
  }


  kingdom_taxa_abundance <-  data.frame()
  unique_kingdom <- unique(kingdom_per_sample$kingdom)
  samples_kingdom <- colnames(kingdom_sample_percent[,4:length(kingdom_sample_percent)])
  for(i in 1:length(samples_kingdom)){
    for(j in 1:length(unique_kingdom)){
      kingdom_taxa_abundance[j,i] <- kingdom_sample_percent %>%
        filter(., kingdom == unique_kingdom[j]) %>%
        select(., samples_kingdom[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(kingdom_taxa_abundance) <- colnames(kingdom_sample_percent[,4:length(kingdom_sample_percent)])

  kingdom_abundance_per_sample <- data.frame(kingdom = unique_kingdom, kingdom_taxa_abundance)
  write_csv(kingdom_abundance_per_sample, file = './kingdom_abundance_per_sample.csv')
  kingdom_abundance_per_sample <- kingdom_abundance_per_sample %>%
    mutate(., kingdom = replace_na(kingdom, 'Rare'))

  return(kingdom_abundance_per_sample)

}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For ASV 18S, using phylum as the taxonomic level.
#' @param asv_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_asv_18s <- function(asv_per_sample){
  asv_full_abundance <-  asv_per_sample[,4] %>%
    map_dbl(sum)

  asv_sample_percent <- data.frame(asv = asv_per_sample$otu,
                                   size = asv_per_sample$size,
                                   genus = asv_per_sample$phylum)
  for(i in 4:length(asv_per_sample)){
    asv_sample_percent[i] <- (asv_per_sample[i]/asv_full_abundance)*100
  }


  asv_taxa_abundance <-  data.frame()
  unique_asv <- unique(asv_per_sample$otu)
  samples_asv <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])
  for(i in 1:length(samples_asv)){
    for(j in 1:length(unique_asv)){
      asv_taxa_abundance[j,i] <- asv_sample_percent %>%
        filter(., asv == unique_asv[j]) %>%
        select(., samples_asv[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(asv_taxa_abundance) <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])

  asv_abundance_per_sample <- data.frame(asv = unique_asv, phylum=asv_per_sample$phylum,
                                         asv_taxa_abundance)
  write_csv(asv_abundance_per_sample, file = './asv_abundance_per_sample.csv')
  asv_abundance_per_sample <- asv_abundance_per_sample %>%
    mutate(., asv = replace_na(asv, 'Rare'))

  return(asv_abundance_per_sample)

}

##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance for each taxa in each sample.
#' NOTE: only works for rarefied data.
#' The relative abundance is calculated based on the total count of the sample on column . Therefore, the calculation is not relative to each sample's total count.
#' For ASV 16S, using genus as the taxonomic level.
#' @param asv_per_sample,  taxonomy cleaned (3-columned object) joined with your metadata
#' @export
get_rel_abun_per_sample_asv_16s <- function(asv_per_sample){
  asv_full_abundance <-  asv_per_sample[,4] %>%
    map_dbl(sum)

  asv_sample_percent <- data.frame(asv = asv_per_sample$otu,
                                   size = asv_per_sample$size,
                                   genus = asv_per_sample$genus)
  for(i in 4:length(asv_per_sample)){
    asv_sample_percent[i] <- (asv_per_sample[i]/asv_full_abundance)*100
  }


  asv_taxa_abundance <-  data.frame()
  unique_asv <- unique(asv_per_sample$otu)
  samples_asv <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])
  for(i in 1:length(samples_asv)){
    for(j in 1:length(unique_asv)){
      asv_taxa_abundance[j,i] <- asv_sample_percent %>%
        filter(., asv == unique_asv[j]) %>%
        select(., samples_asv[i]) %>%
        map_dbl(sum)

    }
  }
  colnames(asv_taxa_abundance) <- colnames(asv_sample_percent[,4:length(asv_sample_percent)])

  asv_abundance_per_sample <- data.frame(asv = unique_asv, genus=asv_per_sample$genus,
                                         asv_taxa_abundance)
  write_csv(asv_abundance_per_sample, file = './asv_abundance_per_sample.csv')
  asv_abundance_per_sample <- asv_abundance_per_sample %>%
    mutate(., asv = replace_na(asv, 'Rare'))

  return(asv_abundance_per_sample)

}

#GET ABUNDANCE TABLE FOR GRAPHING####

##domain####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For domain
#' @param samples_domain, Sample names.
#' @param unique_domain, Each taxon found at this taxonomic level.
#' @param domain_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_domain_abund_graph <-  function(samples_domain, unique_domain,
                                    domain_abundance_per_sample){
  domain_abundance_graph <- data.frame(matrix(NA,nrow = (length(unique_domain)*length(samples_domain)),
                                              ncol = 3))
  #this is for domain binding
  domain_loop <- c()
  for(i in 1:length(unique_domain)){
    for(j in 1:length(samples_domain)){
      domain_loop[i+j] <- (i*j)+1
    }
  }
  domain_loop[length(unique_domain)] <- 1
  domain_loop <- domain_loop[length(unique_domain):(length(domain_loop)-1)]

  for(i in domain_loop){
    for(j in 1:length(unique_domain)){
      domain_abundance_graph[(i+j-1),1] <- unique_domain[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in domain_loop){
    counter <-  counter +1
    for(j in 1:length(unique_domain)){
      domain_abundance_graph[(i+j-1),2] <- samples_domain[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in domain_loop){
    counter <-  counter +1
    for(j in 1:length(unique_domain)){
      domain_abundance_graph[(i+j-1),3] <- domain_abundance_per_sample[j, counter]
    }
  }

  colnames(domain_abundance_graph) <- c('domain', 'sample', 'relative_abundance')

  return(domain_abundance_graph)

}


##class####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For class.
#' @param samples_class, Sample names.
#' @param unique_class, Each taxon found at this taxonomic level.
#' @param class_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_class_abund_graph <-  function(samples_class, unique_class,
                                   class_abundance_per_sample){
  class_abundance_graph <- data.frame(matrix(NA,
                                             nrow = (length(unique_class)*length(samples_class)),
                                             ncol = 3))
  #this is for class binding
  class_loop <- c()
  for(i in 1:length(unique_class)){
    for(j in 1:length(samples_class)){
      class_loop[i+j] <- (i*j)+1
    }
  }
  class_loop[length(unique_class)] <- 1
  class_loop <- class_loop[length(unique_class):(length(class_loop)-1)]

  for(i in class_loop){
    for(j in 1:length(unique_class)){
      class_abundance_graph[(i+j-1),1] <- unique_class[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in class_loop){
    counter <-  counter +1
    for(j in 1:length(unique_class)){
      class_abundance_graph[(i+j-1),2] <- samples_class[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in class_loop){
    counter <-  counter +1
    for(j in 1:length(unique_class)){
      class_abundance_graph[(i+j-1),3] <- class_abundance_per_sample[j, counter]
    }
  }
  colnames(class_abundance_graph) <- c('class', 'sample', 'relative_abundance')

  return(class_abundance_graph)

}


##supergroup####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For supergroup
#' @param samples_supergroup, Sample names.
#' @param unique_supergroup, Each taxon found at this taxonomic level.
#' @param supergroup_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_supergroup_abund_graph <-  function(samples_supergroup, unique_supergroup,
                                        supergroup_abundance_per_sample){
  supergroup_abundance_graph <- data.frame(matrix(NA,nrow = (length(unique_supergroup)*length(samples_supergroup)),
                                                  ncol = 3))
  #this is for supergroup binding
  supergroup_loop <- c()
  for(i in 1:length(unique_supergroup)){
    for(j in 1:length(samples_supergroup)){
      supergroup_loop[i+j] <- (i*j)+1
    }
  }
  supergroup_loop[length(unique_supergroup)] <- 1
  supergroup_loop <- supergroup_loop[length(unique_supergroup):(length(supergroup_loop)-1)]

  for(i in supergroup_loop){
    for(j in 1:length(unique_supergroup)){
      supergroup_abundance_graph[(i+j-1),1] <- unique_supergroup[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in supergroup_loop){
    counter <-  counter +1
    for(j in 1:length(unique_supergroup)){
      supergroup_abundance_graph[(i+j-1),2] <- samples_supergroup[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in supergroup_loop){
    counter <-  counter +1
    for(j in 1:length(unique_supergroup)){
      supergroup_abundance_graph[(i+j-1),3] <- supergroup_abundance_per_sample[j, counter]
    }
  }

  colnames(supergroup_abundance_graph) <- c('supergroup', 'sample', 'relative_abundance')

  return(supergroup_abundance_graph)

}


##division####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For division.
#' @param samples_division, Sample names.
#' @param unique_division, Each taxon found at this taxonomic level.
#' @param division_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_division_abund_graph <-  function(samples_division, unique_division,
                                      division_abundance_per_sample){
  division_abundance_graph <- data.frame(matrix(NA,nrow = (length(unique_division)*length(samples_division)),
                                                ncol = 3))
  #this is for division binding
  division_loop <- c()
  for(i in 1:length(unique_division)){
    for(j in 1:length(samples_division)){
      division_loop[i+j] <- (i*j)+1
    }
  }
  division_loop[length(unique_division)] <- 1
  division_loop <- division_loop[length(unique_division):(length(division_loop)-1)]

  for(i in division_loop){
    for(j in 1:length(unique_division)){
      division_abundance_graph[(i+j-1),1] <- unique_division[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in division_loop){
    counter <-  counter +1
    for(j in 1:length(unique_division)){
      division_abundance_graph[(i+j-1),2] <- samples_division[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in division_loop){
    counter <-  counter +1
    for(j in 1:length(unique_division)){
      division_abundance_graph[(i+j-1),3] <- division_abundance_per_sample[j, counter]
    }
  }

  colnames(division_abundance_graph) <- c('division', 'sample', 'relative_abundance')

  return(division_abundance_graph)

}


##subdivision####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For subdivision.
#' @param samples_subdivision, Sample names.
#' @param unique_subdivision, Each taxon found at this taxonomic level.
#' @param subdivision_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_subdivision_abund_graph <-  function(samples_subdivision, unique_subdivision,
                                         subdivision_abundance_per_sample){
  subdivision_abundance_graph <- data.frame(matrix(NA,nrow = (length(unique_subdivision)*length(samples_subdivision)),
                                                   ncol = 3))
  #this is for subdivision binding
  subdivision_loop <- c()
  for(i in 1:length(unique_subdivision)){
    for(j in 1:length(samples_subdivision)){
      subdivision_loop[i+j] <- (i*j)+1
    }
  }
  subdivision_loop[length(unique_subdivision)] <- 1
  subdivision_loop <- subdivision_loop[length(unique_subdivision):(length(subdivision_loop)-1)]

  for(i in subdivision_loop){
    for(j in 1:length(unique_subdivision)){
      subdivision_abundance_graph[(i+j-1),1] <- unique_subdivision[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in subdivision_loop){
    counter <-  counter +1
    for(j in 1:length(unique_subdivision)){
      subdivision_abundance_graph[(i+j-1),2] <- samples_subdivision[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in subdivision_loop){
    counter <-  counter +1
    for(j in 1:length(unique_subdivision)){
      subdivision_abundance_graph[(i+j-1),3] <- subdivision_abundance_per_sample[j, counter]
    }
  }

  colnames(subdivision_abundance_graph) <- c('subdivision', 'sample', 'relative_abundance')

  return(subdivision_abundance_graph)

}

##order####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For order.
#' @param samples_order, Sample names.
#' @param unique_order, Each taxon found at this taxonomic level.
#' @param order_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_order_abund_graph <-  function(samples_order, unique_order,
                                   order_abundance_per_sample){
  order_abundance_graph <- data.frame(matrix(NA,nrow = (length(unique_order)*length(samples_order)),
                                             ncol = 3))
  #this is for order binding
  order_loop <- c()
  for(i in 1:length(unique_order)){
    for(j in 1:length(samples_order)){
      order_loop[i+j] <- (i*j)+1
    }
  }
  order_loop[length(unique_order)] <- 1
  order_loop <- order_loop[length(unique_order):(length(order_loop)-1)]

  for(i in order_loop){
    for(j in 1:length(unique_order)){
      order_abundance_graph[(i+j-1),1] <- unique_order[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in order_loop){
    counter <-  counter +1
    for(j in 1:length(unique_order)){
      order_abundance_graph[(i+j-1),2] <- samples_order[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in order_loop){
    counter <-  counter +1
    for(j in 1:length(unique_order)){
      order_abundance_graph[(i+j-1),3] <- order_abundance_per_sample[j, counter]
    }
  }

  colnames(order_abundance_graph) <- c('order', 'sample', 'relative_abundance')

  return(order_abundance_graph)

}

##family####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For family.
#' @param samples_family, Sample names.
#' @param unique_family, Each taxon found at this taxonomic level.
#' @param family_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_family_abund_graph <-  function(samples_family, unique_family,
                                    family_abundance_per_sample){
  family_abundance_graph <- data.frame(matrix(NA,
                                              nrow = (length(unique_family)*length(samples_family)),
                                              ncol = 3))
  #this is for family binding
  family_loop <- c()
  for(i in 1:length(unique_family)){
    for(j in 1:length(samples_family)){
      family_loop[i+j] <- (i*j)+1
    }
  }
  family_loop[length(unique_family)] <- 1
  family_loop <- family_loop[length(unique_family):(length(family_loop)-1)]

  for(i in family_loop){
    for(j in 1:length(unique_family)){
      family_abundance_graph[(i+j-1),1] <- unique_family[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in family_loop){
    counter <-  counter +1
    for(j in 1:length(unique_family)){
      family_abundance_graph[(i+j-1),2] <- samples_family[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in family_loop){
    counter <-  counter +1
    for(j in 1:length(unique_family)){
      family_abundance_graph[(i+j-1),3] <- family_abundance_per_sample[j, counter]
    }
  }

  colnames(family_abundance_graph) <- c('family', 'sample', 'relative_abundance')

  return(family_abundance_graph)

}

##genus####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For genus
#' @param samples_genus, Sample names.
#' @param unique_genus, Each taxon found at this taxonomic level.
#' @param genus_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_genus_abund_graph <-  function(samples_genus, unique_genus,
                                   genus_abundance_per_sample){
  genus_abundance_graph <- data.frame(matrix(NA,
                                             nrow = (length(unique_genus)*length(samples_genus)),
                                             ncol = 3))
  #this is for genus binding
  genus_loop <- c()
  for(i in 1:length(unique_genus)){
    for(j in 1:length(samples_genus)){
      genus_loop[i+j] <- (i*j)+1
    }
  }
  genus_loop[length(unique_genus)] <- 1
  genus_loop <- genus_loop[length(unique_genus):(length(genus_loop)-1)]

  for(i in genus_loop){
    for(j in 1:length(unique_genus)){
      genus_abundance_graph[(i+j-1),1] <- unique_genus[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in genus_loop){
    counter <-  counter +1
    for(j in 1:length(unique_genus)){
      genus_abundance_graph[(i+j-1),2] <- samples_genus[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in genus_loop){
    counter <-  counter +1
    for(j in 1:length(unique_genus)){
      genus_abundance_graph[(i+j-1),3] <- genus_abundance_per_sample[j, counter]
    }
  }

  colnames(genus_abundance_graph) <- c('genus', 'sample', 'relative_abundance')

  return(genus_abundance_graph)

}

##species####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For species
#' @param samples_species, Sample names.
#' @param unique_species, Each taxon found at this taxonomic level.
#' @param species_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_species_abund_graph <-  function(samples_species, unique_species,
                                   species_abundance_per_sample){
  species_abundance_graph <- data.frame(matrix(NA,
                                             nrow = (length(unique_species)*length(samples_species)),
                                             ncol = 3))
  #this is for species binding
  species_loop <- c()
  for(i in 1:length(unique_species)){
    for(j in 1:length(samples_species)){
      species_loop[i+j] <- (i*j)+1
    }
  }
  species_loop[length(unique_species)] <- 1
  species_loop <- species_loop[length(unique_species):(length(species_loop)-1)]

  for(i in species_loop){
    for(j in 1:length(unique_species)){
      species_abundance_graph[(i+j-1),1] <- unique_species[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in species_loop){
    counter <-  counter +1
    for(j in 1:length(unique_species)){
      species_abundance_graph[(i+j-1),2] <- samples_species[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in species_loop){
    counter <-  counter +1
    for(j in 1:length(unique_species)){
      species_abundance_graph[(i+j-1),3] <- species_abundance_per_sample[j, counter]
    }
  }

  colnames(species_abundance_graph) <- c('species', 'sample', 'relative_abundance')

  return(species_abundance_graph)

}


##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For ASV 16S, suing genus as a taxonomic level.
#' @param samples_asv, Sample names.
#' @param unique_asv, Each taxon found at this taxonomic level.
#' @param asv_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_asv_abund_graph_16s <-  function(samples_asv, unique_asv,
                                     asv_abundance_per_sample){
  asv_abundance_graph <- data.frame(matrix(NA,
                                           nrow = (length(unique_asv)*length(samples_asv)),
                                           ncol = 4))
  #this is for asv binding
  asv_loop <- c()
  for(i in 1:length(unique_asv)){
    for(j in 1:length(samples_asv)){
      asv_loop[i+j] <- (i*j)+1
    }
  }
  asv_loop[length(unique_asv)] <- 1
  asv_loop <- asv_loop[length(unique_asv):(length(asv_loop)-1)]

  for(i in asv_loop){
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),1] <- unique_asv[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),2] <- samples_asv[counter]
    }
  }

  #this is to add the genus
  counter <-1
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),3] <- asv_abundance_per_sample[j, 2]
    }
  }

  #this is to add the abundance data
  counter <-2
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),4] <- asv_abundance_per_sample[j, counter]
    }
  }

  colnames(asv_abundance_graph) <- c('asv', 'sample', 'genus', 'relative_abundance')

  return(asv_abundance_graph)

}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For ASV 18S, using phylum as the taxonomic level.
#' @param samples_asv, Sample names.
#' @param unique_asv, Each taxon found at this taxonomic level.
#' @param asv_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_asv_abund_graph_18s <-  function(samples_asv, unique_asv,
                                     asv_abundance_per_sample){
  asv_abundance_graph <- data.frame(matrix(NA,
                                           nrow = (length(unique_asv)*length(samples_asv)),
                                           ncol = 4))
  #this is for asv binding
  asv_loop <- c()
  for(i in 1:length(unique_asv)){
    for(j in 1:length(samples_asv)){
      asv_loop[i+j] <- (i*j)+1
    }
  }
  asv_loop[length(unique_asv)] <- 1
  asv_loop <- asv_loop[length(unique_asv):(length(asv_loop)-1)]

  for(i in asv_loop){
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),1] <- unique_asv[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),2] <- samples_asv[counter]
    }
  }

  #this is to add the phylum
  counter <-1
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),3] <- asv_abundance_per_sample[j, 2]
    }
  }

  #this is to add the abundance data
  counter <-2
  for(i in asv_loop){
    counter <-  counter +1
    for(j in 1:length(unique_asv)){
      asv_abundance_graph[(i+j-1),4] <- asv_abundance_per_sample[j, counter]
    }
  }

  colnames(asv_abundance_graph) <- c('asv', 'sample', 'phylum', 'relative_abundance')

  return(asv_abundance_graph)

}

##species####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For generic taxonomic level.
#' @param samples_taxon, Sample names.
#' @param unique_taxon, Each taxon found at this taxonomic level.
#' @param taxon_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_taxon_abund_graph <-  function(samples_taxon, unique_taxon,
                                     taxon_abundance_per_sample){
  taxon_abundance_graph <- data.frame(matrix(NA,
                                               nrow = (length(unique_taxon)*length(samples_taxon)),
                                               ncol = 3))
  #this is for taxon binding
  taxon_loop <- c()
  for(i in 1:length(unique_taxon)){
    for(j in 1:length(samples_taxon)){
      taxon_loop[i+j] <- (i*j)+1
    }
  }
  taxon_loop[length(unique_taxon)] <- 1
  taxon_loop <- taxon_loop[length(unique_taxon):(length(taxon_loop)-1)]

  for(i in taxon_loop){
    for(j in 1:length(unique_taxon)){
      taxon_abundance_graph[(i+j-1),1] <- unique_taxon[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in taxon_loop){
    counter <-  counter +1
    for(j in 1:length(unique_taxon)){
      taxon_abundance_graph[(i+j-1),2] <- samples_taxon[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in taxon_loop){
    counter <-  counter +1
    for(j in 1:length(unique_taxon)){
      taxon_abundance_graph[(i+j-1),3] <- taxon_abundance_per_sample[j, counter]
    }
  }

  colnames(taxon_abundance_graph) <- c('taxon', 'sample', 'relative_abundance')

  return(taxon_abundance_graph)

}

##kingdom####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For Kingdom.
#' @param samples_kingdom, Sample names.
#' @param unique_kingdom, Each taxon found at this taxonomic level.
#' @param kingdom_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_kingdom_abund_graph <-  function(samples_kingdom, unique_kingdom,
                                     kingdom_abundance_per_sample){
  kingdom_abundance_graph <- data.frame(matrix(NA,
                                               nrow = (length(unique_kingdom)*length(samples_kingdom)),
                                               ncol = 3))
  #this is for kingdom binding
  kingdom_loop <- c()
  for(i in 1:length(unique_kingdom)){
    for(j in 1:length(samples_kingdom)){
      kingdom_loop[i+j] <- (i*j)+1
    }
  }
  kingdom_loop[length(unique_kingdom)] <- 1
  kingdom_loop <- kingdom_loop[length(unique_kingdom):(length(kingdom_loop)-1)]

  for(i in kingdom_loop){
    for(j in 1:length(unique_kingdom)){
      kingdom_abundance_graph[(i+j-1),1] <- unique_kingdom[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in kingdom_loop){
    counter <-  counter +1
    for(j in 1:length(unique_kingdom)){
      kingdom_abundance_graph[(i+j-1),2] <- samples_kingdom[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in kingdom_loop){
    counter <-  counter +1
    for(j in 1:length(unique_kingdom)){
      kingdom_abundance_graph[(i+j-1),3] <- kingdom_abundance_per_sample[j, counter]
    }
  }

  colnames(kingdom_abundance_graph) <- c('kingdom', 'sample', 'relative_abundance')

  return(kingdom_abundance_graph)

}


##phylum####
#' @title Microinformatics #' @description
#' This function is to reshape the *taxa*_abundance_per_sample object, in which each column is a sample,
#' to just 3 columns: taxa, sample, abundance.
#' For Phylum.
#' @param samples_phylum, Sample names.
#' @param unique_phylum, Each taxon found at this taxonomic level.
#' @param phylum_abundance_per_sample, Table where each taxon is a row and each column is a sample containing the total or relative abundance
#' @export
get_phylum_abund_graph <-  function(samples_phylum, unique_phylum,
                                    phylum_abundance_per_sample){
  phylum_abundance_graph <- data.frame(matrix(NA,
                                              nrow = (length(unique_phylum)*length(samples_phylum)),
                                              ncol = 3))
  #this is for phylum binding
  phylum_loop <- c()
  for(i in 1:length(unique_phylum)){
    for(j in 1:length(samples_phylum)){
      phylum_loop[i+j] <- (i*j)+1
    }
  }
  phylum_loop[length(unique_phylum)] <- 1
  phylum_loop <- phylum_loop[length(unique_phylum):(length(phylum_loop)-1)]

  for(i in phylum_loop){
    for(j in 1:length(unique_phylum)){
      phylum_abundance_graph[(i+j-1),1] <- unique_phylum[j]
    }
  }


  #this is for sample binding


  counter <-0
  for(i in phylum_loop){
    counter <-  counter +1
    for(j in 1:length(unique_phylum)){
      phylum_abundance_graph[(i+j-1),2] <- samples_phylum[counter]
    }
  }

  #this is to add the abundance data
  counter <-1
  for(i in phylum_loop){
    counter <-  counter +1
    for(j in 1:length(unique_phylum)){
      phylum_abundance_graph[(i+j-1),3] <- phylum_abundance_per_sample[j, counter]
    }
  }

  colnames(phylum_abundance_graph) <- c('phylum', 'sample', 'relative_abundance')

  return(phylum_abundance_graph)

}

#get taxa abundance for graphing####

##domain####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For domain.
#' @param filter_domain, Taxonomy table cleaned in a 3-column format (ASV, size, domain)
#' @export
get_domain_taxa_abundance <- function(filter_domain){
  unique_domain <- unique(filter_domain$domain)
  domain_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_domain),
                                              ncol = 1))
  for(j in 1:length(unique_domain)){
    domain_taxa_abundance[j,1] <- filter_domain %>%
      filter(., domain == unique_domain[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(domain_taxa_abundance) <- 'size'
  domain_taxa_abundance <- data_frame(domain = unique_domain,
                                      domain_taxa_abundance)
  return(domain_taxa_abundance)
}


##supergroup####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For supergroup.
#' @param filter_supergroup, Taxonomy table cleaned in a 3-column format (ASV, size, supergroup)
#' @export
get_supergroup_taxa_abundance <- function(filter_supergroup){
  unique_supergroup <- unique(filter_supergroup$supergroup)
  supergroup_taxa_abundance <-  data.frame(matrix(NA,
                                                  nrow = length(unique_supergroup),
                                                  ncol = 1))
  for(j in 1:length(unique_supergroup)){
    supergroup_taxa_abundance[j,1] <- filter_supergroup %>%
      filter(., supergroup == unique_supergroup[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(supergroup_taxa_abundance) <- 'size'
  supergroup_taxa_abundance <- data_frame(supergroup = unique_supergroup,
                                          supergroup_taxa_abundance)
  return(supergroup_taxa_abundance)
}

##division####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For division.
#' @param filter_division, Taxonomy table cleaned in a 3-column format (ASV, size, division)
#' @export
get_division_taxa_abundance <- function(filter_division){
  unique_division <- unique(filter_division$division)
  division_taxa_abundance <-  data.frame(matrix(NA,
                                                nrow = length(unique_division),
                                                ncol = 1))
  for(j in 1:length(unique_division)){
    division_taxa_abundance[j,1] <- filter_division %>%
      filter(., division == unique_division[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(division_taxa_abundance) <- 'size'
  division_taxa_abundance <- data_frame(division = unique_division,
                                        division_taxa_abundance)
  return(division_taxa_abundance)
}

##subdivision####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For subdivision.
#' @param filter_subdivision, Taxonomy table cleaned in a 3-column format (ASV, size, subdivision)
#' @export
get_subdivision_taxa_abundance <- function(filter_subdivision){
  unique_subdivision <- unique(filter_subdivision$subdivision)
  subdivision_taxa_abundance <-  data.frame(matrix(NA,
                                                   nrow = length(unique_subdivision),
                                                   ncol = 1))
  for(j in 1:length(unique_subdivision)){
    subdivision_taxa_abundance[j,1] <- filter_subdivision %>%
      filter(., subdivision == unique_subdivision[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(subdivision_taxa_abundance) <- 'size'
  subdivision_taxa_abundance <- data_frame(subdivision = unique_subdivision,
                                           subdivision_taxa_abundance)
  return(subdivision_taxa_abundance)
}

##class####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For class.
#' @param filter_class, Taxonomy table cleaned in a 3-column format (ASV, size, class)
#' @export

get_class_taxa_abundance <- function(filter_class){
  unique_class <- unique(filter_class$class)
  class_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_class),
                                             ncol = 1))
  for(j in 1:length(unique_class)){
    class_taxa_abundance[j,1] <- filter_class %>%
      filter(., class == unique_class[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(class_taxa_abundance) <- 'size'
  class_taxa_abundance <- data_frame(class = unique_class,
                                     class_taxa_abundance)
  return(class_taxa_abundance)
}

##order####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For order.
#' @param filter_order, Taxonomy table cleaned in a 3-column format (ASV, size, order)
#' @export
get_order_taxa_abundance <- function(filter_order){
  unique_order <- unique(filter_order$order)
  order_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_order),
                                             ncol = 1))
  for(j in 1:length(unique_order)){
    order_taxa_abundance[j,1] <- filter_order %>%
      filter(., order == unique_order[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(order_taxa_abundance) <- 'size'
  order_taxa_abundance <- data_frame(order = unique_order,
                                     order_taxa_abundance)
  return(order_taxa_abundance)
}


##family####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For family.
#' @param filter_family, Taxonomy table cleaned in a 3-column format (ASV, size, family)
#' @export
get_family_taxa_abundance <- function(filter_family){
  unique_family <- unique(filter_family$family)
  family_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_family),
                                              ncol = 1))
  for(j in 1:length(unique_family)){
    family_taxa_abundance[j,1] <- filter_family %>%
      filter(., family == unique_family[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(family_taxa_abundance) <- 'size'
  family_taxa_abundance <- data_frame(family = unique_family,
                                      family_taxa_abundance)
  return(family_taxa_abundance)
}

##genus####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For genus.
#' @param filter_genus, Taxonomy table cleaned in a 3-column format (ASV, size, genus)
#' @export
get_genus_taxa_abundance <- function(filter_genus){
  unique_genus <- unique(filter_genus$genus)
  genus_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_genus),
                                             ncol = 1))
  for(j in 1:length(unique_genus)){
    genus_taxa_abundance[j,1] <- filter_genus %>%
      filter(., genus == unique_genus[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(genus_taxa_abundance) <- 'size'
  genus_taxa_abundance <- data_frame(genus = unique_genus,
                                     genus_taxa_abundance)
  return(genus_taxa_abundance)
}

##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For ASV 16S, using genus as the taxonomic level.
#' @param filter_genus, Taxonomy table cleaned in a 3-column format (ASV, size, genus)
#' @export
get_asv_taxa_abundance_16s <- function(filter_genus){
  unique_asv <- unique(filter_genus$otu)
  asv_taxa_abundance <-  data.frame(matrix(NA,
                                           nrow = length(unique_asv),
                                           ncol = 1))
  for(j in 1:length(unique_asv)){
    asv_taxa_abundance[j,1] <- filter_genus %>%
      filter(., otu == unique_asv[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(asv_taxa_abundance) <- 'size'
  asv_taxa_abundance <- data_frame(asv = unique_asv, genus = filter_genus$genus,
                                   asv_taxa_abundance)
  return(asv_taxa_abundance)
}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For ASV 18S, using phylum as the taxonomic reference level.
#' @param filter_phylum, Taxonomy table cleaned in a 3-column format (ASV, size, phylum)
#' @export
get_asv_taxa_abundance_18s <- function(filter_phylum){
  unique_asv <- unique(filter_phylum$otu)
  asv_taxa_abundance <-  data.frame(matrix(NA,
                                           nrow = length(unique_asv),
                                           ncol = 1))
  for(j in 1:length(unique_asv)){
    asv_taxa_abundance[j,1] <- filter_phylum %>%
      filter(., otu == unique_asv[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(asv_taxa_abundance) <- 'size'
  asv_taxa_abundance <- data_frame(asv = unique_asv, phylum = filter_phylum$phylum,
                                   asv_taxa_abundance)
  return(asv_taxa_abundance)
}

##phylum####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For phylum.
#' @param filter_phylum, Taxonomy table cleaned in a 3-column format (ASV, size, phylum)
#' @export
get_phylum_taxa_abundance <- function(filter_phylum){
  unique_phylum <- unique(filter_phylum$phylum)
  phylum_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_phylum),
                                              ncol = 1))
  for(j in 1:length(unique_phylum)){
    phylum_taxa_abundance[j,1] <- filter_phylum %>%
      filter(., phylum == unique_phylum[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(phylum_taxa_abundance) <- 'size'
  phylum_taxa_abundance <- data_frame(phylum = unique_phylum,
                                      phylum_taxa_abundance)
  return(phylum_taxa_abundance)
}

##kingdom####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For kingdom.
#' @param filter_kingdom, Taxonomy table cleaned in a 3-column format (ASV, size, kingdom)
#' @export
get_kingdom_taxa_abundance <- function(filter_kingdom){
  unique_kingdom <- unique(filter_kingdom$kingdom)
  kingdom_taxa_abundance <-  data.frame(matrix(NA,
                                               nrow = length(unique_kingdom),
                                               ncol = 1))
  for(j in 1:length(unique_kingdom)){
    kingdom_taxa_abundance[j,1] <- filter_kingdom %>%
      filter(., kingdom == unique_kingdom[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(kingdom_taxa_abundance) <- 'size'
  kingdom_taxa_abundance <- data_frame(kingdom = unique_kingdom,
                                       kingdom_taxa_abundance)
  return(kingdom_taxa_abundance)
}

##species####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For species.
#' @param filter_species, Taxonomy table cleaned in a 3-column format (ASV, size, species)
#' @export
get_species_taxa_abundance <- function(filter_species){
  unique_species <- unique(filter_species$species)
  species_taxa_abundance <-  data.frame(matrix(NA,
                                               nrow = length(unique_species),
                                               ncol = 1))
  for(j in 1:length(unique_species)){
    species_taxa_abundance[j,1] <- filter_species %>%
      filter(., species == unique_species[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(species_taxa_abundance) <- 'size'
  species_taxa_abundance <- data_frame(species = unique_species,
                                       species_taxa_abundance)
  return(species_taxa_abundance)
}

##taxon####
#' @title Microinformatics #' @description
#' This function is to get the total abundance for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For non-specified taxonomic level
#' @param filter_taxon, Taxonomy table cleaned in a 3-column format (ASV, size, taxon)
#' @export
get_taxon_taxa_abundance <- function(filter_taxon){
  unique_taxon <- unique(pull(filter_taxon[,3]))
  taxon_taxa_abundance <-  data.frame(matrix(NA,
                                               nrow = length(unique_taxon),
                                               ncol = 1))
  for(j in 1:length(unique_taxon)){
    taxon_taxa_abundance[j,1] <- filter_taxon %>%
      filter(., taxon == unique_taxon[j]) %>%
      select(., size) %>%
      map_dbl(sum)
  }
  colnames(taxon_taxa_abundance) <- 'size'
  taxon_taxa_abundance <- data_frame(taxon = unique_taxon,
                                       taxon_taxa_abundance)
  return(taxon_taxa_abundance)
}


#get_taxa_abundance_per_sample####

##domain####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For domain.
#' @param domain_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_domain_taxa_abundance_per_sample <- function(domain_per_sample){
  unique_domain <- unique(domain_per_sample$domain)
  domain_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_domain),
                                              ncol = ncol(domain_per_sample)-3))
  for(i in 1:(ncol(domain_per_sample)-3)){
    domain_per_sample_counter <- domain_per_sample %>%
      select(., otu, size, domain, i+3)
    for(j in 1:length(unique_domain)){
      domain_taxa_abundance[j,i] <- domain_per_sample_counter %>%
        filter(., domain == unique_domain[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(domain_taxa_abundance) <- colnames(domain_per_sample[,4:ncol(domain_per_sample)])
  domain_taxa_abundance <- data_frame(domain = unique_domain,
                                      domain_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(domain_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(domain_taxa_abundance[i,2:ncol(domain_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  domain_taxa_abundance <- data_frame(domain_taxa_abundance,  total_sum_samples)

  return(domain_taxa_abundance)
}


##subdivision####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For subdivision.
#' @param subdivision_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_subdivision_taxa_abundance_per_sample <- function(subdivision_per_sample){
  unique_subdivision <- unique(subdivision_per_sample$subdivision)
  subdivision_taxa_abundance <-  data.frame(matrix(NA,
                                                   nrow = length(unique_subdivision),
                                                   ncol = ncol(subdivision_per_sample)-3))
  for(i in 1:(ncol(subdivision_per_sample)-3)){
    subdivision_per_sample_counter <- subdivision_per_sample %>%
      select(., otu, size, subdivision, i+3)
    for(j in 1:length(unique_subdivision)){
      subdivision_taxa_abundance[j,i] <- subdivision_per_sample_counter %>%
        filter(., subdivision == unique_subdivision[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(subdivision_taxa_abundance) <- colnames(subdivision_per_sample[,4:ncol(subdivision_per_sample)])
  subdivision_taxa_abundance <- data_frame(subdivision = unique_subdivision,
                                           subdivision_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(subdivision_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(subdivision_taxa_abundance[i,2:ncol(subdivision_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  subdivision_taxa_abundance <- data_frame(subdivision_taxa_abundance,  total_sum_samples)

  return(subdivision_taxa_abundance)
}


##division####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For division.
#' @param division_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_division_taxa_abundance_per_sample <- function(division_per_sample){
  unique_division <- unique(division_per_sample$division)
  division_taxa_abundance <-  data.frame(matrix(NA,
                                                nrow = length(unique_division),
                                                ncol = ncol(division_per_sample)-3))
  for(i in 1:(ncol(division_per_sample)-3)){
    division_per_sample_counter <- division_per_sample %>%
      select(., otu, size, division, i+3)
    for(j in 1:length(unique_division)){
      division_taxa_abundance[j,i] <- division_per_sample_counter %>%
        filter(., division == unique_division[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(division_taxa_abundance) <- colnames(division_per_sample[,4:ncol(division_per_sample)])
  division_taxa_abundance <- data_frame(division = unique_division,
                                        division_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(division_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(division_taxa_abundance[i,2:ncol(division_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  division_taxa_abundance <- data_frame(division_taxa_abundance,  total_sum_samples)

  return(division_taxa_abundance)
}

##supergroup####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For supergroup.
#' @param supergroup_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_supergroup_taxa_abundance_per_sample <- function(supergroup_per_sample){
  unique_supergroup <- unique(supergroup_per_sample$supergroup)
  supergroup_taxa_abundance <-  data.frame(matrix(NA,
                                                  nrow = length(unique_supergroup),
                                                  ncol = ncol(supergroup_per_sample)-3))
  for(i in 1:(ncol(supergroup_per_sample)-3)){
    supergroup_per_sample_counter <- supergroup_per_sample %>%
      select(., otu, size, supergroup, i+3)
    for(j in 1:length(unique_supergroup)){
      supergroup_taxa_abundance[j,i] <- supergroup_per_sample_counter %>%
        filter(., supergroup == unique_supergroup[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(supergroup_taxa_abundance) <- colnames(supergroup_per_sample[,4:ncol(supergroup_per_sample)])
  supergroup_taxa_abundance <- data_frame(supergroup = unique_supergroup,
                                          supergroup_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(supergroup_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(supergroup_taxa_abundance[i,2:ncol(supergroup_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  supergroup_taxa_abundance <- data_frame(supergroup_taxa_abundance,  total_sum_samples)

  return(supergroup_taxa_abundance)
}

##phylum####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For phylum.
#' @param phylum_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_phylum_taxa_abundance_per_sample <- function(phylum_per_sample){
  unique_phylum <- unique(phylum_per_sample$phylum)
  phylum_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_phylum),
                                              ncol = ncol(phylum_per_sample)-3))
  for(i in 1:(ncol(phylum_per_sample)-3)){
    phylum_per_sample_counter <- phylum_per_sample %>%
      select(., otu, size, phylum, i+3)
    for(j in 1:length(unique_phylum)){
      phylum_taxa_abundance[j,i] <- phylum_per_sample_counter %>%
        filter(., phylum == unique_phylum[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(phylum_taxa_abundance) <- colnames(phylum_per_sample[,4:ncol(phylum_per_sample)])
  phylum_taxa_abundance <- data_frame(phylum = unique_phylum,
                                      phylum_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(phylum_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(phylum_taxa_abundance[i,2:ncol(phylum_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  phylum_taxa_abundance <- data_frame(phylum_taxa_abundance,  total_sum_samples)

  return(phylum_taxa_abundance)
}


##class####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For class.
#' @param class_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_class_taxa_abundance_per_sample <- function(class_per_sample){
  unique_class <- unique(class_per_sample$class)
  class_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_class),
                                             ncol = ncol(class_per_sample)-3))
  for(i in 1:(ncol(class_per_sample)-3)){
    class_per_sample_counter <- class_per_sample %>%
      select(., otu, size, class, i+3)
    for(j in 1:length(unique_class)){
      class_taxa_abundance[j,i] <- class_per_sample_counter %>%
        filter(., class == unique_class[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(class_taxa_abundance) <- colnames(class_per_sample[,4:ncol(class_per_sample)])
  class_taxa_abundance <- data_frame(class = unique_class,
                                     class_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(class_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(class_taxa_abundance[i,2:ncol(class_taxa_abundance)])
  }
  colnames(total_sum_samples) <- c('total_sum')
  class_taxa_abundance <- data_frame(class_taxa_abundance, total_sum_samples)

  return(class_taxa_abundance)
}

##order####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For order.
#' @param order_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_order_taxa_abundance_per_sample <- function(order_per_sample){
  unique_order <- unique(order_per_sample$order)
  order_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_order),
                                             ncol = ncol(order_per_sample)-3))
  for(i in 1:(ncol(order_per_sample)-3)){
    order_per_sample_counter <- order_per_sample %>%
      select(., otu, size, order, i+3)
    for(j in 1:length(unique_order)){
      order_taxa_abundance[j,i] <- order_per_sample_counter %>%
        filter(., order == unique_order[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(order_taxa_abundance) <- colnames(order_per_sample[,4:ncol(order_per_sample)])
  order_taxa_abundance <- data_frame(order = unique_order,
                                     order_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(order_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(order_taxa_abundance[i,2:ncol(order_taxa_abundance)])
  }
  colnames(total_sum_samples) <- c('total_sum')
  order_taxa_abundance <- data_frame(order_taxa_abundance, total_sum_samples)

  return(order_taxa_abundance)
}

##family####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For family.
#' @param family_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_family_taxa_abundance_per_sample <- function(family_per_sample){
  unique_family <- unique(family_per_sample$family)
  family_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_family),
                                              ncol = ncol(family_per_sample)-3))
  for(i in 1:(ncol(family_per_sample)-3)){
    family_per_sample_counter <- family_per_sample %>%
      select(., otu, size, family, i+3)
    for(j in 1:length(unique_family)){
      family_taxa_abundance[j,i] <- family_per_sample_counter %>%
        filter(., family == unique_family[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(family_taxa_abundance) <- colnames(family_per_sample[,4:ncol(family_per_sample)])
  family_taxa_abundance <- data_frame(family = unique_family,
                                      family_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(family_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(family_taxa_abundance[i,2:ncol(family_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')

  family_taxa_abundance <- data_frame(family_taxa_abundance, total_sum_samples)

  return(family_taxa_abundance)
}

##genus####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For genus.
#' @param genus_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_genus_taxa_abundance_per_sample <- function(genus_per_sample){
  unique_genus <- unique(genus_per_sample$genus)
  genus_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_genus),
                                             ncol = ncol(genus_per_sample)-3))
  for(i in 1:(ncol(genus_per_sample)-3)){
    genus_per_sample_counter <- genus_per_sample %>%
      select(., otu, size, genus, i+3)
    for(j in 1:length(unique_genus)){
      genus_taxa_abundance[j,i] <- genus_per_sample_counter %>%
        filter(., genus == unique_genus[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(genus_taxa_abundance) <- colnames(genus_per_sample[,4:ncol(genus_per_sample)])
  genus_taxa_abundance <- data_frame(genus = unique_genus,
                                     genus_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(genus_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(genus_taxa_abundance[i,2:ncol(genus_taxa_abundance)])
  }
  colnames(total_sum_samples) <- c('total_sum')
  genus_taxa_abundance <- data_frame(genus_taxa_abundance, total_sum_samples)

  return(genus_taxa_abundance)
}

##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For ASV 16S, using genus as the reference taxonomic level.
#' @param genus_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_asv_16s_taxa_abundance_per_sample <- function(genus_per_sample){
  unique_asv <- unique(genus_per_sample$otu)
  genus_taxa_abundance <-  data.frame(matrix(NA,
                                             nrow = length(unique_asv),
                                             ncol = ncol(genus_per_sample)-3))
  for(i in 1:(ncol(genus_per_sample)-3)){
    genus_per_sample_counter <- genus_per_sample %>%
      select(., otu, size, genus, i+3)
    for(j in 1:length(unique_asv)){
      genus_taxa_abundance[j,i] <- genus_per_sample_counter %>%
        filter(., otu == unique_asv[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(genus_taxa_abundance) <- colnames(genus_per_sample[,4:ncol(genus_per_sample)])
  genus_taxa_abundance <- data_frame(asv = unique_asv, genus = genus_per_sample$genus,
                                     genus_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(genus_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(genus_taxa_abundance[i,3:ncol(genus_taxa_abundance)])
  }
  colnames(total_sum_samples) <- c('total_sum')
  genus_taxa_abundance <- data_frame(genus_taxa_abundance, total_sum_samples)

  return(genus_taxa_abundance)
}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For ASV 18S, using phylum as the reference taxonomic level.
#' @param phylum_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_asv_18s_taxa_abundance_per_sample <- function(phylum_per_sample){
  unique_asv <- unique(phylum_per_sample$otu)
  phylum_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_asv),
                                              ncol = ncol(phylum_per_sample)-3))
  for(i in 1:(ncol(phylum_per_sample)-3)){
    phylum_per_sample_counter <- phylum_per_sample %>%
      select(., otu, size, phylum, i+3)
    for(j in 1:length(unique_asv)){
      phylum_taxa_abundance[j,i] <- phylum_per_sample_counter %>%
        filter(., otu == unique_asv[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(phylum_taxa_abundance) <- colnames(phylum_per_sample[,4:ncol(phylum_per_sample)])
  phylum_taxa_abundance <- data_frame(asv = unique_asv, phylum = phylum_per_sample$phylum,
                                      phylum_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(phylum_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(phylum_taxa_abundance[i,3:ncol(phylum_taxa_abundance)])
  }
  colnames(total_sum_samples) <- c('total_sum')
  phylum_taxa_abundance <- data_frame(phylum_taxa_abundance, total_sum_samples)

  return(phylum_taxa_abundance)
}

##species####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For species.
#' @param species_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_species_taxa_abundance_per_sample <- function(species_per_sample){
  unique_species <- unique(species_per_sample$species)
  species_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_species),
                                              ncol = ncol(species_per_sample)-3))
  for(i in 1:(ncol(species_per_sample)-3)){
    species_per_sample_counter <- species_per_sample %>%
      select(., otu, size, species, i+3)
    for(j in 1:length(unique_species)){
      species_taxa_abundance[j,i] <- species_per_sample_counter %>%
        filter(., species == unique_species[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(species_taxa_abundance) <- colnames(species_per_sample[,4:ncol(species_per_sample)])
  species_taxa_abundance <- data_frame(species = unique_species,
                                      species_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(species_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(species_taxa_abundance[i,2:ncol(species_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  species_taxa_abundance <- data_frame(species_taxa_abundance,  total_sum_samples)

  return(species_taxa_abundance)
}

##taxon####
#' @title Microinformatics #' @description
#' This function is to get the total abundance of each taxa for each sample.
#' It also adds a final column with the total count sum for each taxa across all samples.
#' The function is done for each taxonomic level.
#' For non-specific taxonomic level.
#' @param taxon_per_sample, Taxonomy cleaned joined with the shared table having samples as columns.
#' @export
get_taxon_taxa_abundance_per_sample <- function(taxon_per_sample){
  unique_taxon <- unique(pull(taxon_per_sample[,3]))
  taxon_taxa_abundance <-  data.frame(matrix(NA,
                                              nrow = length(unique_taxon),
                                              ncol = ncol(taxon_per_sample)-3))
  for(i in 1:(ncol(taxon_per_sample)-3)){
    taxon_per_sample_counter <- taxon_per_sample %>%
      select(., otu, size, taxon, i+3)
    for(j in 1:length(unique_taxon)){
      taxon_taxa_abundance[j,i] <- taxon_per_sample_counter %>%
        filter(., taxon == unique_taxon[j]) %>%
        select(., 4) %>%
        map_dbl(sum)
    }
  }


  colnames(taxon_taxa_abundance) <- colnames(taxon_per_sample[,4:ncol(taxon_per_sample)])
  taxon_taxa_abundance <- data_frame(taxon = unique_taxon,
                                      taxon_taxa_abundance)
  total_sum_samples <- data.frame()
  for(i in 1:nrow(taxon_taxa_abundance)){
    total_sum_samples[i,1] <-  rowSums(taxon_taxa_abundance[i,2:ncol(taxon_taxa_abundance)])
  }

  colnames(total_sum_samples) <- c('total_sum')
  taxon_taxa_abundance <- data_frame(taxon_taxa_abundance,  total_sum_samples)

  return(taxon_taxa_abundance)
}



#GET TAXA OTHERS RELATIVE ABUNDANCE####

##domain####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For domain.
#' @param domain_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_domain, The samples' names.
#' @param top_taxa_domain, The most abundant taxon by abundance previously defined by the user.
#' @export
get_domain_all_others <- function(domain_for_graph, samples_domain,
                                  top_taxa_domain){
  domain_graph_all_others <- tibble()
  domain_looper <- tibble()
  for(i in 1:length(samples_domain)){
    domain_looper <- domain_for_graph %>%
      select(., domain, sample, relative_abundance) %>%
      filter(., sample == samples_domain[i]) %>%
      filter(., !(domain %in% top_taxa_domain))
    domain_graph_all_others <- rbind(domain_graph_all_others, domain_looper)

  }
  return(domain_graph_all_others)
}

##supergroup####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For supergroup
#' @param supergroup_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_supergroup, The samples' names.
#' @param top_taxa_supergroup, The most abundant taxon by abundance previously defined by the user.
#' @export
get_supergroup_all_others <- function(supergroup_for_graph, samples_supergroup,
                                      top_taxa_supergroup){
  supergroup_graph_all_others <- tibble()
  supergroup_looper <- tibble()
  for(i in 1:length(samples_supergroup)){
    supergroup_looper <- supergroup_for_graph %>%
      select(., supergroup, sample, relative_abundance) %>%
      filter(., sample == samples_supergroup[i]) %>%
      filter(., !(supergroup %in% top_taxa_supergroup))
    supergroup_graph_all_others <- rbind(supergroup_graph_all_others, supergroup_looper)

  }
  return(supergroup_graph_all_others)
}

##division####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For division.
#' @param division_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_division, The samples' names.
#' @param top_taxa_division, The most abundant taxon by abundance previously defined by the user.
#' @export
get_division_all_others <- function(division_for_graph, samples_division,
                                    top_taxa_division){
  division_graph_all_others <- tibble()
  division_looper <- tibble()
  for(i in 1:length(samples_division)){
    division_looper <- division_for_graph %>%
      select(., division, sample, relative_abundance) %>%
      filter(., sample == samples_division[i]) %>%
      filter(., !(division %in% top_taxa_division))
    division_graph_all_others <- rbind(division_graph_all_others, division_looper)

  }
  return(division_graph_all_others)
}

##subdivision####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For subdivision.
#' @param subdivision_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_subdivision, The samples' names.
#' @param top_taxa_subdivision, The most abundant taxon by abundance previously defined by the user.
#' @export
get_subdivision_all_others <- function(subdivision_for_graph, samples_subdivision,
                                       top_taxa_subdivision){
  subdivision_graph_all_others <- tibble()
  subdivision_looper <- tibble()
  for(i in 1:length(samples_subdivision)){
    subdivision_looper <- subdivision_for_graph %>%
      select(., subdivision, sample, relative_abundance) %>%
      filter(., sample == samples_subdivision[i]) %>%
      filter(., !(subdivision %in% top_taxa_subdivision))
    subdivision_graph_all_others <- rbind(subdivision_graph_all_others, subdivision_looper)

  }
  return(subdivision_graph_all_others)
}

##kingdom####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For kingdom.
#' @param kingdom_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_kingdom, The samples' names.
#' @param top_taxa_kingdom, The most abundant taxon by abundance previously defined by the user.
#' @export
get_kingdom_all_others <- function(kingdom_for_graph, samples_kingdom,
                                   top_taxa_kingdom){
  kingdom_graph_all_others <- tibble()
  kingdom_looper <- tibble()
  for(i in 1:length(samples_kingdom)){
    kingdom_looper <- kingdom_for_graph %>%
      select(., kingdom, sample, relative_abundance) %>%
      filter(., sample == samples_kingdom[i]) %>%
      filter(., !(kingdom %in% top_taxa_kingdom))
    kingdom_graph_all_others <- rbind(kingdom_graph_all_others, kingdom_looper)

  }
  return(kingdom_graph_all_others)
}

##phylum####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For phylum.
#' @param phylum_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_phylum, The samples' names.
#' @param top_taxa_phylum, The most abundant taxon by abundance previously defined by the user.
#' @export
get_phylum_all_others <- function(phylum_for_graph, samples_phylum,
                                  top_taxa_phylum){
  phylum_graph_all_others <- tibble()
  phylum_looper <- tibble()
  for(i in 1:length(samples_phylum)){
    phylum_looper <- phylum_for_graph %>%
      select(., phylum, sample, relative_abundance) %>%
      filter(., sample == samples_phylum[i]) %>%
      filter(., !(phylum %in% top_taxa_phylum))
    phylum_graph_all_others <- rbind(phylum_graph_all_others, phylum_looper)

  }
  return(phylum_graph_all_others)
}

##class####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For class.
#' @param class_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_class, The samples' names.
#' @param top_taxa_class, The most abundant taxon by abundance previously defined by the user.
#' @export
get_class_all_others <- function(class_for_graph, samples_class,
                                 top_taxa_class){
  class_graph_all_others <- tibble()
  class_looper <- tibble()
  for(i in 1:length(samples_class)){
    class_looper <- class_for_graph %>%
      select(., class, sample, relative_abundance) %>%
      filter(., sample == samples_class[i]) %>%
      filter(., !(class %in% top_taxa_class))
    class_graph_all_others <- rbind(class_graph_all_others, class_looper)

  }
  return(class_graph_all_others)
}

##order####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For order.
#' @param order_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_order, The samples' names.
#' @param top_taxa_order, The most abundant taxon by abundance previously defined by the user.
#' @export
get_order_all_others <- function(order_for_graph, samples_order,
                                 top_taxa_order){
  order_graph_all_others <- tibble()
  order_looper <- tibble()
  for(i in 1:length(samples_order)){
    order_looper <- order_for_graph %>%
      select(., order, sample, relative_abundance) %>%
      filter(., sample == samples_order[i]) %>%
      filter(., !(order %in% top_taxa_order))
    order_graph_all_others <- rbind(order_graph_all_others, order_looper)

  }
  return(order_graph_all_others)
}

##family####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For family.
#' @param family_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_family, The samples' names.
#' @param top_taxa_family, The most abundant taxon by abundance previously defined by the user.
#' @export
get_family_all_others <- function(family_for_graph, samples_family,
                                  top_taxa_family){
  family_graph_all_others <- tibble()
  family_looper <- tibble()
  for(i in 1:length(samples_family)){
    family_looper <- family_for_graph %>%
      select(., family, sample, relative_abundance) %>%
      filter(., sample == samples_family[i]) %>%
      filter(., !(family %in% top_taxa_family))
    family_graph_all_others <- rbind(family_graph_all_others, family_looper)

  }
  return(family_graph_all_others)
}

##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For ASV 16s, using genus as the reference taxonomic level.
#' @param asv_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_asv, The samples' names.
#' @param top_taxa_asv, The most abundant taxon by abundance previously defined by the user.
#' @export
get_asv_all_others_16s <- function(asv_for_graph, samples_asv,
                                   top_taxa_asv){
  asv_graph_all_others <- tibble()
  asv_looper <- tibble()
  for(i in 1:length(samples_asv)){
    asv_looper <- asv_for_graph %>%
      select(., asv, genus, sample, relative_abundance) %>%
      filter(., sample == samples_asv[i]) %>%
      filter(., !(asv %in% top_taxa_asv))
    asv_graph_all_others <- rbind(asv_graph_all_others, asv_looper)

  }
  return(asv_graph_all_others)
}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For ASV 18s, using phylum as the reference taxonomic level.
#' @param asv_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_asv, The samples' names.
#' @param top_taxa_asv, The most abundant taxon by abundance previously defined by the user.
#' @export
get_asv_all_others_18s <- function(asv_for_graph, samples_asv,
                                   top_taxa_asv){
  asv_graph_all_others <- tibble()
  asv_looper <- tibble()
  for(i in 1:length(samples_asv)){
    asv_looper <- asv_for_graph %>%
      select(., asv, phylum, sample, relative_abundance) %>%
      filter(., sample == samples_asv[i]) %>%
      filter(., !(asv %in% top_taxa_asv))
    asv_graph_all_others <- rbind(asv_graph_all_others, asv_looper)

  }
  return(asv_graph_all_others)
}

##genus####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For genus
#' @param genus_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_genus, The samples' names.
#' @param top_taxa_genus, The most abundant taxon by abundance previously defined by the user.
#' @export
get_genus_all_others <- function(genus_for_graph, samples_genus,
                                 top_taxa_genus){
  genus_graph_all_others <- tibble()
  genus_looper <- tibble()
  for(i in 1:length(samples_genus)){
    genus_looper <- genus_for_graph %>%
      select(., genus, sample, relative_abundance) %>%
      filter(., sample == samples_genus[i]) %>%
      filter(., !(genus %in% top_taxa_genus))
    genus_graph_all_others <- rbind(genus_graph_all_others, genus_looper)

  }
  return(genus_graph_all_others)
}

##taxon####
#' @title Microinformatics #' @description
#' This function is to get the relative abundance of taxa that are not part of the Top X ranking.
#' The user must have previously determine the most abundance taxa (i.e. the top 15 taxa).
#' The function retains 3 columns with the taxa, sample and relative abundance of all the other taxa.
#' For non-specified taxonomic level.
#' @param taxon_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_taxon, The samples' names.
#' @param top_taxa_taxon, The most abundant taxon by abundance previously defined by the user.
#' @export
get_taxon_all_others <- function(taxon_for_graph, samples_taxon,
                                 top_taxa_taxon){
  taxon_graph_all_others <- tibble()
  taxon_looper <- tibble()
  for(i in 1:length(samples_taxon)){
    taxon_looper <- taxon_for_graph %>%
      select(., taxon, sample, relative_abundance) %>%
      filter(., sample == samples_taxon[i]) %>%
      filter(., !(taxon %in% top_taxa_taxon))
    taxon_graph_all_others <- rbind(taxon_graph_all_others, taxon_looper)

  }
  return(taxon_graph_all_others)
}


#get taxa for graph with others####
##genus####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For genus.
#' @param genus_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_genus, The samples' names.
#' @param top_taxa_genus, The most abundant taxon by abundance previously defined by the user.
#' @export
get_genus_sum_others <- function(genus_for_graph, samples_genus,
                                 top_taxa_genus){
  genus_graph_sum_others <- tibble()
  genus_looper <- tibble()
  genus_looper_3 <- tibble(sample=genus_for_graph[,2], genus_for_graph[,4:ncol(genus_for_graph)])
  for(i in 1:length(samples_genus)){
    genus_looper <- genus_for_graph %>%
      select(., genus, sample, relative_abundance) %>%
      filter(., sample == samples_genus[i]) %>%
      filter(., !(genus %in% top_taxa_genus)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    genus_looper_2 <- tibble(genus = c('Others'), sample = c(samples_genus[i]),
                             relative_abundance = c(genus_looper))
    genus_looper_4 <-inner_join(genus_looper_2, genus_looper_3,
                                by = c('sample'='sample'))
    genus_looper_4 <- genus_looper_4[1,1:51]
    genus_graph_sum_others <- rbind(genus_graph_sum_others, genus_looper_4)
  }
  return(genus_graph_sum_others)
}

##taxon####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For non-specified taxonomic level.
#' @param taxon_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_taxon, The samples' names.
#' @param top_taxa_taxon, The most abundant taxon by abundance previously defined by the user.
#' @export
get_taxon_sum_others <- function(taxon_for_graph, samples_taxon,
                                 top_taxa_taxon){
  taxon_graph_sum_others <- tibble()
  taxon_looper <- tibble()
  taxon_looper_3 <- tibble(sample=taxon_for_graph[,2], genus_for_graph[,4:ncol(taxon_for_graph)])
  for(i in 1:length(samples_taxon)){
    taxon_looper <- taxon_for_graph %>%
      select(., taxon, sample, relative_abundance) %>%
      filter(., sample == samples_taxon[i]) %>%
      filter(., !(taxon %in% top_taxa_taxon)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    taxon_looper_2 <- tibble(taxon = c('Others'), sample = c(samples_taxon[i]),
                             relative_abundance = c(taxon_looper))
    taxon_looper_4 <-inner_join(taxon_looper_2, taxon_looper_3,
                                by = c('sample'='sample'))
    taxon_looper_4 <- taxon_looper_4[1,1:51]
    taxon_graph_sum_others <- rbind(taxon_graph_sum_others, taxon_looper_4)
  }
  return(taxon_graph_sum_others)
}

##phylum####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For phylum.
#' @param phylum_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_phylum, The samples' names.
#' @param top_taxa_phylum, The most abundant taxon by abundance previously defined by the user.
#' @export
get_phylum_sum_others <- function(phylum_for_graph, samples_phylum,
                                  top_taxa_phylum){
  phylum_graph_sum_others <- tibble()
  phylum_looper <- tibble()
  phylum_looper_3 <- tibble(sample=phylum_for_graph[,2], phylum_for_graph[,4:ncol(phylum_for_graph)])
  for(i in 1:length(samples_phylum)){
    phylum_looper <- phylum_for_graph %>%
      select(., phylum, sample, relative_abundance) %>%
      filter(., sample == samples_phylum[i]) %>%
      filter(., !(phylum %in% top_taxa_phylum)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    phylum_looper_2 <- tibble(phylum = c('Others'), sample = c(samples_phylum[i]),
                              relative_abundance = c(phylum_looper))
    phylum_looper_4 <-inner_join(phylum_looper_2, phylum_looper_3,
                                 by = c('sample'='sample'))
    phylum_looper_4 <- phylum_looper_4[1,1:ncol(phylum_for_graph)]
    phylum_graph_sum_others <- rbind(phylum_graph_sum_others, phylum_looper_4)
  }
  return(phylum_graph_sum_others)
}

##order####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For order.
#' @param order_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_order, The samples' names.
#' @param top_taxa_order, The most abundant taxon by abundance previously defined by the user.
#' @export
get_order_sum_others <- function(order_for_graph, samples_order,
                                  top_taxa_order){
  order_graph_sum_others <- tibble()
  order_looper <- tibble()
  order_looper_3 <- tibble(sample=order_for_graph[,2], order_for_graph[,4:ncol(order_for_graph)])
  for(i in 1:length(samples_order)){
    order_looper <- order_for_graph %>%
      select(., order, sample, relative_abundance) %>%
      filter(., sample == samples_order[i]) %>%
      filter(., !(order %in% top_taxa_order)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    order_looper_2 <- tibble(order = c('Others'), sample = c(samples_order[i]),
                              relative_abundance = c(order_looper))
    order_looper_4 <-inner_join(order_looper_2, order_looper_3,
                                 by = c('sample'='sample'))
    order_looper_4 <- order_looper_4[1,1:ncol(order_for_graph)]
    order_graph_sum_others <- rbind(order_graph_sum_others, order_looper_4)
  }
  return(order_graph_sum_others)
}

##class####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For class.
#' @param class_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_class, The samples' names.
#' @param top_taxa_class, The most abundant taxon by abundance previously defined by the user.
#' @export
get_class_sum_others <- function(class_for_graph, samples_class,
                                  top_taxa_class){
  class_graph_sum_others <- tibble()
  class_looper <- tibble()
  class_looper_3 <- tibble(sample=class_for_graph[,2], class_for_graph[,4:ncol(class_for_graph)])
  for(i in 1:length(samples_class)){
    class_looper <- class_for_graph %>%
      select(., class, sample, relative_abundance) %>%
      filter(., sample == samples_class[i]) %>%
      filter(., !(class %in% top_taxa_class)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    class_looper_2 <- tibble(class = c('Others'), sample = c(samples_class[i]),
                              relative_abundance = c(class_looper))
    class_looper_4 <-inner_join(class_looper_2, class_looper_3,
                                 by = c('sample'='sample'))
    class_looper_4 <- class_looper_4[1,1:ncol(class_for_graph)]
    class_graph_sum_others <- rbind(class_graph_sum_others, class_looper_4)
  }
  return(class_graph_sum_others)
}

##family####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For family.
#' @param family_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_family, The samples' names.
#' @param top_taxa_family, The most abundant taxon by abundance previously defined by the user.
#' @export
get_family_sum_others <- function(family_for_graph, samples_family,
                                  top_taxa_family){
  family_graph_sum_others <- tibble()
  family_looper <- tibble()
  family_looper_3 <- tibble(sample=family_for_graph[,2], family_for_graph[,4:ncol(family_for_graph)])
  for(i in 1:length(samples_family)){
    family_looper <- family_for_graph %>%
      select(., family, sample, relative_abundance) %>%
      filter(., sample == samples_family[i]) %>%
      filter(., !(family %in% top_taxa_family)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    family_looper_2 <- tibble(family = c('Others'), sample = c(samples_family[i]),
                              relative_abundance = c(family_looper))
    family_looper_4 <-inner_join(family_looper_2, family_looper_3,
                                 by = c('sample'='sample'))
    family_looper_4 <- family_looper_4[1,1:ncol(family_for_graph)]
    family_graph_sum_others <- rbind(family_graph_sum_others, family_looper_4)
  }
  return(family_graph_sum_others)
}

##species####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For species.
#' @param species_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_species, The samples' names.
#' @param top_taxa_species, The most abundant taxon by abundance previously defined by the user.
#' @export
get_species_sum_others <- function(species_for_graph, samples_species,
                                  top_taxa_species){
  species_graph_sum_others <- tibble()
  species_looper <- tibble()
  species_looper_3 <- tibble(sample=species_for_graph[,2], species_for_graph[,4:ncol(species_for_graph)])
  for(i in 1:length(samples_species)){
    species_looper <- species_for_graph %>%
      select(., species, sample, relative_abundance) %>%
      filter(., sample == samples_species[i]) %>%
      filter(., !(species %in% top_taxa_species)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    species_looper_2 <- tibble(species = c('Others'), sample = c(samples_species[i]),
                              relative_abundance = c(species_looper))
    species_looper_4 <-inner_join(species_looper_2, species_looper_3,
                                 by = c('sample'='sample'))
    species_looper_4 <- species_looper_4[1,1:ncol(species_for_graph)]
    species_graph_sum_others <- rbind(species_graph_sum_others, species_looper_4)
  }
  return(species_graph_sum_others)
}

##kingdom####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For kingdom.
#' @param kingdom_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_kingdom, The samples' names.
#' @param top_taxa_kingdom, The most abundant taxon by abundance previously defined by the user.
#' @export
get_kingdom_sum_others <- function(kingdom_for_graph, samples_kingdom,
                                  top_taxa_kingdom){
  kingdom_graph_sum_others <- tibble()
  kingdom_looper <- tibble()
  kingdom_looper_3 <- tibble(sample=kingdom_for_graph[,2], kingdom_for_graph[,4:ncol(kingdom_for_graph)])
  for(i in 1:length(samples_kingdom)){
    kingdom_looper <- kingdom_for_graph %>%
      select(., kingdom, sample, relative_abundance) %>%
      filter(., sample == samples_kingdom[i]) %>%
      filter(., !(kingdom %in% top_taxa_kingdom)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    kingdom_looper_2 <- tibble(kingdom = c('Others'), sample = c(samples_kingdom[i]),
                              relative_abundance = c(kingdom_looper))
    kingdom_looper_4 <-inner_join(kingdom_looper_2, kingdom_looper_3,
                                 by = c('sample'='sample'))
    kingdom_looper_4 <- kingdom_looper_4[1,1:ncol(kingdom_for_graph)]
    kingdom_graph_sum_others <- rbind(kingdom_graph_sum_others, kingdom_looper_4)
  }
  return(kingdom_graph_sum_others)
}

##domain####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For domain.
#' @param domain_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_domain, The samples' names.
#' @param top_taxa_domain, The most abundant taxon by abundance previously defined by the user.
#' @export
get_domain_sum_others <- function(domain_for_graph, samples_domain,
                                  top_taxa_domain){
  domain_graph_sum_others <- tibble()
  domain_looper <- tibble()
  domain_looper_3 <- tibble(sample=domain_for_graph[,2], domain_for_graph[,4:ncol(domain_for_graph)])
  for(i in 1:length(samples_domain)){
    domain_looper <- domain_for_graph %>%
      select(., domain, sample, relative_abundance) %>%
      filter(., sample == samples_domain[i]) %>%
      filter(., !(domain %in% top_taxa_domain)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    domain_looper_2 <- tibble(domain = c('Others'), sample = c(samples_domain[i]),
                              relative_abundance = c(domain_looper))
    domain_looper_4 <-inner_join(domain_looper_2, domain_looper_3,
                                 by = c('sample'='sample'))
    domain_looper_4 <- domain_looper_4[1,1:ncol(domain_for_graph)]
    domain_graph_sum_others <- rbind(domain_graph_sum_others, domain_looper_4)
  }
  return(domain_graph_sum_others)
}

##supergroup####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For supergroup.
#' @param supergroup_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_supergroup, The samples' names.
#' @param top_taxa_supergroup, The most abundant taxon by abundance previously defined by the user.
#' @export
get_supergroup_sum_others <- function(supergroup_for_graph, samples_supergroup,
                                  top_taxa_supergroup){
  supergroup_graph_sum_others <- tibble()
  supergroup_looper <- tibble()
  supergroup_looper_3 <- tibble(sample=supergroup_for_graph[,2], supergroup_for_graph[,4:ncol(supergroup_for_graph)])
  for(i in 1:length(samples_supergroup)){
    supergroup_looper <- supergroup_for_graph %>%
      select(., supergroup, sample, relative_abundance) %>%
      filter(., sample == samples_supergroup[i]) %>%
      filter(., !(supergroup %in% top_taxa_supergroup)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    supergroup_looper_2 <- tibble(supergroup = c('Others'), sample = c(samples_supergroup[i]),
                              relative_abundance = c(supergroup_looper))
    supergroup_looper_4 <-inner_join(supergroup_looper_2, supergroup_looper_3,
                                 by = c('sample'='sample'))
    supergroup_looper_4 <- supergroup_looper_4[1,1:ncol(supergroup_for_graph)]
    supergroup_graph_sum_others <- rbind(supergroup_graph_sum_others, supergroup_looper_4)
  }
  return(supergroup_graph_sum_others)
}

##division####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For division.
#' @param division_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_division, The samples' names.
#' @param top_taxa_division, The most abundant taxon by abundance previously defined by the user.
#' @export
get_division_sum_others <- function(division_for_graph, samples_division,
                                  top_taxa_division){
  division_graph_sum_others <- tibble()
  division_looper <- tibble()
  division_looper_3 <- tibble(sample=division_for_graph[,2], division_for_graph[,4:ncol(division_for_graph)])
  for(i in 1:length(samples_division)){
    division_looper <- division_for_graph %>%
      select(., division, sample, relative_abundance) %>%
      filter(., sample == samples_division[i]) %>%
      filter(., !(division %in% top_taxa_division)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    division_looper_2 <- tibble(division = c('Others'), sample = c(samples_division[i]),
                              relative_abundance = c(division_looper))
    division_looper_4 <-inner_join(division_looper_2, division_looper_3,
                                 by = c('sample'='sample'))
    division_looper_4 <- division_looper_4[1,1:ncol(division_for_graph)]
    division_graph_sum_others <- rbind(division_graph_sum_others, division_looper_4)
  }
  return(division_graph_sum_others)
}

##subdivision####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For subdivision.
#' @param subdivision_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_subdivision, The samples' names.
#' @param top_taxa_subdivision, The most abundant taxon by abundance previously defined by the user.
#' @export
get_subdivision_sum_others <- function(subdivision_for_graph, samples_subdivision,
                                  top_taxa_subdivision){
  subdivision_graph_sum_others <- tibble()
  subdivision_looper <- tibble()
  subdivision_looper_3 <- tibble(sample=subdivision_for_graph[,2], subdivision_for_graph[,4:ncol(subdivision_for_graph)])
  for(i in 1:length(samples_subdivision)){
    subdivision_looper <- subdivision_for_graph %>%
      select(., subdivision, sample, relative_abundance) %>%
      filter(., sample == samples_subdivision[i]) %>%
      filter(., !(subdivision %in% top_taxa_subdivision)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    subdivision_looper_2 <- tibble(subdivision = c('Others'), sample = c(samples_subdivision[i]),
                              relative_abundance = c(subdivision_looper))
    subdivision_looper_4 <-inner_join(subdivision_looper_2, subdivision_looper_3,
                                 by = c('sample'='sample'))
    subdivision_looper_4 <- subdivision_looper_4[1,1:ncol(subdivision_for_graph)]
    subdivision_graph_sum_others <- rbind(subdivision_graph_sum_others, subdivision_looper_4)
  }
  return(subdivision_graph_sum_others)
}

##ASV 16S####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For ASV 16S, using genus as the reference taxonomic level.
#' @param asv_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_asv, The samples' names.
#' @param top_taxa_asv, The most abundant taxon by abundance previously defined by the user.
#' @export
get_asv_sum_others_16s <- function(asv_for_graph, samples_asv,
                                   top_taxa_asv){
  asv_graph_sum_others <- tibble()
  asv_looper <- tibble()
  asv_looper_3 <- tibble(sample=asv_for_graph[,2], genus = asv_for_graph[,3],
                         asv_for_graph[,5:43])
  for(i in 1:length(samples_asv)){
    asv_looper <- asv_for_graph %>%
      select(., asv, genus, sample, relative_abundance) %>%
      filter(., sample == samples_asv[i]) %>%
      filter(., !(asv %in% top_taxa_asv)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    asv_looper_2 <- tibble(asv = c('Others'), sample = c(samples_asv[i]),
                           relative_abundance = c(asv_looper))
    asv_looper_4 <-inner_join(asv_looper_2, asv_looper_3,
                              by = c('sample'='sample'))
    asv_looper_4 <- asv_looper_4[1,1:43]
    asv_graph_sum_others <- rbind(asv_graph_sum_others, asv_looper_4)
  }
  return(asv_graph_sum_others)
}

##ASV 18S####
#' @title Microinformatics #' @description
#' This function is to make an object with the relative abundance, sample, taxa and metadata for those taxa that are not within the "Top abundant" taxa.
#' The top abundant taxa must be previously determine by the user (i.e. top 15 taxa).
#' For ASV 18S, using phylum as the reference taxonomic level.
#' @param asv_for_graph, Table containing the taxonomy cleaned, relative abundance, and metadata.
#' @param samples_asv, The samples' names.
#' @param top_taxa_asv, The most abundant taxon by abundance previously defined by the user.
#' @export
get_asv_sum_others_18s <- function(asv_for_graph, samples_asv,
                                   top_taxa_asv){
  asv_graph_sum_others <- tibble()
  asv_looper <- tibble()
  asv_looper_3 <- tibble(sample=asv_for_graph[,2], phylum = asv_for_graph[,3],
                         asv_for_graph[,5:ncol(asv_for_graph)])
  for(i in 1:length(samples_asv)){
    asv_looper <- asv_for_graph %>%
      select(., asv, phylum, sample, relative_abundance) %>%
      filter(., sample == samples_asv[i]) %>%
      filter(., !(asv %in% top_taxa_asv)) %>%
      select(., relative_abundance) %>%
      map_dbl(sum)
    asv_looper_2 <- tibble(asv = c('Others'), sample = c(samples_asv[i]),
                           relative_abundance = c(asv_looper))
    asv_looper_4 <-inner_join(asv_looper_2, asv_looper_3,
                              by = c('sample'='sample'))
    asv_looper_4 <- asv_looper_4[1,1:ncol(asv_for_graph)]
    asv_graph_sum_others <- rbind(asv_graph_sum_others, asv_looper_4)
  }
  return(asv_graph_sum_others)
}



