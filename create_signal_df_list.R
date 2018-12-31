create_signal_df_list <- function(bin.list, track.list){
  signal.df.list <- list()
  for (i in 1:length(bin.list)){
    gene_bins <- bin.list[[i]]
    gene <- names(bin.list)[i]
    cat(sprintf("\n\n%s\n",gene))
    
    binned_gene_signals <- list()
    for (j in 1:length(track.list)){
      # undebug(bin_signal)
      mark <- names(track.list)[j]
      print(mark)
      # if(gene == "CYP7A1" & mark == "RNAseq_p") debug(bin_signal)
      # binned_gene_signals[[mark]] <- bin_signal(bin.list[[i]], track.list[[j]])
      binned_gene_signals[[mark]] <- mybinnedaverage(bin.list[[i]], track.list[[j]])
    }
    
    #make matrix
    df <- data.frame(mids = start(binned_gene_signals[[1]]) + floor( width(binned_gene_signals[[1]])/2 ))
    for (mark in names(binned_gene_signals)){
      df[,mark] <- as.data.frame(values(binned_gene_signals[[mark]]))
    }
    
    #store in list
    signal.df.list[[gene]] <- df
  }
  return(signal.df.list)
}