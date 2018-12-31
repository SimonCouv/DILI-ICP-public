#generate panel files
#REMARK: panel file already excludes related individuals, eg compare with pedigree file integrated_call_samples.20130502.ALL.ped
  # HG00144 and HG00155 are mother/child (family ID GBR001), only HG00155 in panel file
sample_pop <- fread("../data/1KG_LD/integrated_call_samples_v3.20130502.ALL.panel")
pops_requested <- list(EUR="EUR", GBR="GBR", IBS = "IBS", FIN = "FIN", TSI = "TSI", CEU = "CEU", GBR_CEU = c("GBR", "CEU"))
for (pop_name in names(pops_requested)){
  sample_pop %>% dplyr::filter((super_pop %in% pops_requested[[pop_name]]) | (pop %in% pops_requested[[pop_name]])) %>% 
    dplyr::select(sample) %>%
    # dplyr::select(sample, gender) %>% 
    # dplyr::mutate(gender=dplyr::recode(gender, female = "F", male = "M")) %>%
    write_tsv(.,paste0("../data/1KG_LD/", pop_name, "_samples.tsv"), col_names = FALSE)
}
