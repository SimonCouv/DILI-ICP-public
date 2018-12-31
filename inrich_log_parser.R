
inrich_log_parser <- function(log_fp){
  library(stringr)
  
  log_con <- file(log_fp, open ="r")
  inrich_results <- list()
  while (TRUE){
    # browser()
    outerline <- readLines(log_con, n=1)
    
    if (length(outerline)==0)
      break
    
    if (grepl(outerline, pattern="INRICH")){
      run_block <- c()
      while (TRUE){
        innerline <- readLines(log_con, n=1)
        run_block <- c(run_block,innerline)
        if (grepl(innerline, pattern="Writing output to"))
          break
      }
      
      run_name <- na.omit(str_match(run_block, pattern = "intervals_(\\w{2,3}_\\w{3}_\\d)\\.tsv")[,2])
      cat(sprintf("\nParsing run: %s.\n", run_name))
    
      #enrichment statistics at level of gene set (ICP genes)
      geneset_stat_header_index <- grep(run_block, pattern = "T_Size  Int_No    Empirical_P Corrected_P     Target")
      geneset_stat_line <- run_block[geneset_stat_header_index + 2]
      

    
      inrich_run_df = data_frame(
        dili_pheno = na.omit(str_match(run_block, pattern = "intervals_(\\w{2,3})_\\w{3}_\\d\\.tsv")[,2]),
        population = na.omit(str_match(run_block, pattern = "intervals_\\w{2,3}_(\\w{3})_\\d\\.tsv")[,2]),
        rsq_10 = na.omit(str_match(run_block, pattern = "intervals_\\w{2,3}_\\w{3}_(\\d)\\.tsv")[,2]),
        n_interval_in = na.omit(str_match(run_block, pattern = "read (\\d+) intervals")[,2]),
        n_interval_genic = na.omit(str_match(run_block, pattern = "(\\d+) test elements on gene regions")[,2]),
        n_interval_genic_nonoverlap = na.omit(str_match(run_block, pattern = "After merging, (\\d+) non-overlapping intervals")[,2]),
        
        
        n_overlaps = str_split(geneset_stat_line, pattern = " ")[[1]] %>% .[nchar(.)>0] %>% .[2],
        overlap_p_emp = str_split(geneset_stat_line, pattern = " ")[[1]] %>% .[nchar(.)>0] %>% .[3],
        overlap_p_cor = str_split(geneset_stat_line, pattern = " ")[[1]] %>% .[nchar(.)>0] %>% .[4]
      )
      
      inrich_results <- c(inrich_results, list(inrich_run_df))

    }
    
  }
  inrich_results_df <- bind_rows(inrich_results)
  return(inrich_results_df)
}
