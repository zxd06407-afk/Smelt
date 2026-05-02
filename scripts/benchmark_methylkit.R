#!/usr/bin/env Rscript
# methylKit benchmark for comparison with Smelt
# Arabidopsis met1-3 vs Col-0 WT WGBS data
# methylKit v1.26.0

library(methylKit)

setwd("/home/hermes/Smelt/benchmark_data")
out_dir <- "methylkit_out"
dir.create(out_dir, showWarnings = FALSE)

samples <- c("ERR965674", "ERR965675", "ERR965676", "ERR965677")
# Use CX_report files which contain explicit context information
cx_files <- file.path("smelt_out", paste0(samples, ".CX_report.txt"))

for (ctx in c("CpG", "CHG", "CHH")) {
    cat(sprintf("\n==========================\n=== Context: %s ===\n==========================\n", ctx))

    # 1. Read coverage files
    t1 <- Sys.time()
    meth <- methRead(location = as.list(cx_files),
        sample.id = as.list(samples),
        assembly = "tair10",
        pipeline = "bismarkCytosineReport",
        treatment = c(1, 1, 0, 0),
        context = ctx,
        mincov = 5)
    cat(sprintf("Read: %.0f sec\n", difftime(Sys.time(), t1, units = "secs")))

    # 2. Filter
    meth_f <- filterByCoverage(meth, lo.count = 5, hi.perc = 99.9)

    # 3. Tile + unite
    t2 <- Sys.time()
    tiles <- tileMethylCounts(meth_f, win.size = 2000, step.size = 500)
    meth_unite <- unite(tiles, destrand = TRUE)
    n_win <- if (is(meth_unite, "methylBase")) nrow(meth_unite) else length(meth_unite)
    cat(sprintf("Tiled+unite: %d windows, %.0f sec\n", n_win, difftime(Sys.time(), t2, units = "secs")))

    # 4. DMR
    t3 <- Sys.time()
    dm <- calculateDiffMeth(meth_unite, overdispersion = "none", test = "fast.fisher", mc.cores = 8)
    cat(sprintf("DMR calc: %.0f sec\n", difftime(Sys.time(), t3, units = "secs")))

    dmr <- getMethylDiff(dm, difference = 20, qvalue = 0.05)
    n_dmr <- if (is.null(dmr)) 0 else nrow(dmr)
    cat(sprintf("DMRs q<0.05 |d|>0.2: %d\n", n_dmr))

    # 5. Save
    perc <- percMethylation(meth_unite)
    write.table(perc, file.path(out_dir, sprintf("methylkit_%s_window_meth.tsv", ctx)),
        sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(dm), file.path(out_dir, sprintf("methylkit_%s_dm_all.tsv", ctx)),
        sep = "\t", quote = FALSE, row.names = FALSE)
    if (n_dmr > 0) {
        write.table(as.data.frame(dmr), file.path(out_dir, sprintf("methylkit_%s_dmr.tsv", ctx)),
            sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

cat("\n=== methylKit benchmark complete ===\n")
