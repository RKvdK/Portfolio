# This script retrieves the final output table from the Monte Carlo simulations

MCnum <- 3L
res <- vector("list", nrow(simgrid))

pb <- txtProgressBar(min = 0, max = nrow(simgrid), style = 3)

for (i in seq_len(nrow(simgrid))) {
  g <- simgrid[i, , drop = FALSE]
  
  simout <- SIM(
    N = g$N,
    nA = g$nA,
    nB = g$nB,
    nE = g$nE,
    PercOverlap = g$PercOverlap,
    tranmat_diag = g$tranmat_diag,
    tranmat_sym  = g$tranmat_sym,
    cia = g$cia,
    zpx = g$zpx,
    Xsource = as.character(g$Xsource),
    seed = as.integer(g$ID)
  )
  
  mcout <- MC(SIMout = simout, MCnum = MCnum, seed = as.integer(g$ID))
  
  res[[i]] <- data.frame(
    ID = g$ID,
    mechanism = g$mechanism,
    pattern = g$pattern,
    degree = g$degree,
    DRE_rmse_mean = mcout$DRE$rmse_mean,
    IPF_rmse_mean = mcout$IPF$rmse_mean,
    EXT_rmse_mean = mcout$EXT$rmse_mean,
    MCnum = MCnum
  )
  
  setTxtProgressBar(pb, i)
}

close(pb)

results <- do.call(rbind, res)
write.csv(results, "outputs/grid_results.csv", row.names = FALSE)