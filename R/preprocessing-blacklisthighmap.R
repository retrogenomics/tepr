.wmeanvec <- function(dupframenbvec, currenttrans) {

    wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans) {

        ## Selecting all rows having a window equal to nbdup
        allframedf <- currenttrans[which(currenttrans$window.window == nbdup), ]

        ## Retrieving start and end of the window to calculate nb of nt
        windowstart <- unique(allframedf$start.window)
        windowend <- unique(allframedf$end.window)
        if (!isTRUE(all.equal(length(windowstart), 1)) ||
            !isTRUE(all.equal(length(windowend), 1)))
            stop("\n[tepr] Error: Non-unique window size.\n",
                "  Contact the developer.\n")

        ## Retrieve the nb of overlapping nt for each score
        endselect <- pmin(allframedf$end, windowend)
        startselect <- pmax(allframedf$start, windowstart)
        overntvec <- pmax(endselect - startselect + 1, 0)

        ## Computing weighted mean
        allscores <- as.data.frame(allframedf[, "score"])[[1]]
        wmean <- stats::weighted.mean(allscores, overntvec)
        return(wmean)
    }, currenttrans)

    return(wmeanvec)
}

.removeblackandlowmap <- function(currenttrans, blacklisttib, idxscore,
    maptracktib) {

        ## Set scores overlapping black list to NA
        resblack <- valr::bed_intersect(currenttrans, blacklisttib,
            min_overlap = 0L)
        if (!isTRUE(all.equal(nrow(resblack), 0))) {
            strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                currenttrans$end, sep = "-")
            strblack <- paste(resblack$chrom, resblack$start.x, resblack$end.x,
                sep = "-")
            idxblack <- as.vector(na.omit(unique(match(strtransvec, strblack))))
            if (isTRUE(all.equal(length(idxblack), 0)))
                stop("\n[tepr] Error: Blacklist overlap processing failed.\n",
                    "  Contact the developer.\n")
            currenttrans[idxblack, idxscore] <- NA
        }

        ## Set scores NOT overlapping high map to NA (i.e. scores overlapping
        ## low mappability intervals)
        currenttrans$chrom <- as.character(currenttrans$chrom)
        maptracktib$chrom <- as.character(maptracktib$chrom)
        resmap <- valr::bed_intersect(currenttrans, maptracktib,
            min_overlap = 0L)
        if (!isTRUE(all.equal(nrow(resmap), 0))) {
            ## Compute strtransvec only if no overlap with black list was found
            if (isTRUE(all.equal(nrow(resblack), 0)))
                strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                    currenttrans$end, sep = "-")
            strmap <- paste(resmap$chrom, resmap$start.x, resmap$end.x,
                sep = "-")
            idxmap <- match(strtransvec, strmap)
            idxtosetna <- which(is.na(idxmap)) ## NOT overlapping high map
            currenttrans[idxtosetna, idxscore] <- NA
        }
        return(currenttrans)
}

.meanblackhighbytrans <- function(bgscorebytrans, windsize, currentname,
    currentchrom, blacklisttib, maptracktib, saveobjectpath, nbcputrans,
    reload, verbose) {

        currentobj <- file.path(saveobjectpath, paste0(currentname,
            "-", currentchrom, "-translist.rds"))
        if (!reload || !file.exists(currentobj)) {
            bytranslist <- parallel::mclapply(bgscorebytrans,
                function(currenttrans, windsize, currentname, blacklisttib,
                    maptracktib) {

                        ## Reordering rows according to window number
                        idxwind <- order(currenttrans$window.window)
                        currenttrans <- currenttrans[idxwind, ]

                        ## Identifying duplicated windows that will be used to
                        ## compute a weighted mean.
                        currenttrans <- .dupidx(currenttrans, windsize)

                        ## Remove columns corresponding to bedgraph, move score
                        ## column at the end, remove the .window suffix from
                        ## column names, Add exp name prefix to column score
                        res <- .formatcurrenttranscols(currenttrans,
                            currentname)
                        currenttrans <- res[[1]]
                        idxscore <- res[[2]]

                        ## Set scores overlapping black list and low map to NA
                        idxchrom <- which(maptracktib$chrom == unique(
                            currenttrans$chrom))
                        maptracktibchrom <- maptracktib[idxchrom, ]
                        currenttrans <- .removeblackandlowmap(currenttrans,
                            blacklisttib, idxscore, maptracktibchrom)

                        return(currenttrans)

            }, windsize, currentname, blacklisttib, maptracktib,
                mc.cores = nbcputrans)
            if (!is.na(saveobjectpath)) {
                if (verbose) message("\t\t\t Saving ", currentobj)
                saveRDS(bytranslist, file = currentobj)
            }
        } else {
            if (verbose) message("\t\t\t Loading ", currentobj)
            bytranslist <- readRDS(currentobj)
        }

        return(bytranslist)
}

.retrieveandfilterfrombg <- function(exptab, blacklisttib, maptracktib, # nolint
    nbcputrans, allwindchromtib, expnamevec, windsize, currentchrom,
    chromlength, saveobjectpath, showtime, showmemory, reload, tmpfold,
    verbose) {

        ## Looping on each experiment bg file
        if (verbose) message("\t\t For each bedgraph file") # nolint
        invisible(mapply(function(currentpath, currentname,
            currentstrand, currentcond, currentrep, currentdirection,
            allwindchromtib, blacklisttib, maptracktib, windsize, currentchrom,
            chromlength, nbcputrans, saveobjectpath, verbose, showtime,
            showmemory, reload, tmpfold) {

            filename <- file.path(tmpfold, paste0(currentname, "-",
                currentchrom, ".tsv"))

            if (!reload || !file.exists(filename)) {
                ## Retrieving bedgraph values
                if (verbose) message("\n\t\t Retrieving begraph values for ", # nolint
                    currentname, " on ", currentchrom) # nolint
                valtib <- .retrievebgval(currentpath, currentchrom, chromlength,
                    showmemory, verbose)

                ## Retrieving scores on annotations of strand
                annoscores <- .retrieveannoscores(currentstrand,
                    allwindchromtib, valtib, showmemory, verbose)

                ## Splitting the scores by transcript
                if (verbose) message("\t\t Splitting the scores by transcript")
                trsfact <- factor(annoscores$transcript.window)
                bgscorebytrans <- split(annoscores, trsfact)
                rm(trsfact, annoscores, valtib)
                if (showmemory) print(gc()) else invisible(gc())

                ## For each transcript compute the weighted means for each
                ## window. The weight is calculated if a window contains more
                ## than one score
                if (verbose) message("\t\t For each transcript compute the ",
                    "weighted means and set scores overlapping black list and ",
                    "low mappability to NA. It takes a while.")
                if (showtime) start_time_bytranslist <- Sys.time()

                bytranslist <- .meanblackhighbytrans(bgscorebytrans, windsize,
                    currentname, currentchrom, blacklisttib, maptracktib,
                    saveobjectpath, nbcputrans, reload, verbose)

                if (showtime) {
                    end_time_bytranslist <- Sys.time()
                    timing <- end_time_bytranslist - start_time_bytranslist
                    message("\t\t\t ## Features excluded in: ",
                        format(timing, digits = 2))
                }

                if (!isTRUE(all.equal(unique(sapply(bytranslist, nrow)),
                    windsize)))
                    stop("\n[tepr] Error: Row count mismatch in transcript list.\n",
                        "  Expected ", windsize, " rows per element.\n",
                        "  Contact the developer.\n")

                ## Formatting columns and adding rowid column
                res <- .rowidandcols(bytranslist, currentcond, currentrep,
                    currentdirection, showmemory, verbose)

                ## Saving table to temporary folder
                if (verbose) message("\t\t Saving table to ", filename)
                utils::write.table(res, file = filename, sep = "\t", quote = FALSE,
                    col.names = FALSE, row.names = FALSE)

                rm(bgscorebytrans, bytranslist, res)
                if (showmemory) print(gc()) else invisible(gc())
            } else {
                if (verbose) message("\t\t The file ", filename,
                    " was already computed. Skipping.")
            }

        }, exptab$path, expnamevec, exptab$strand, exptab$condition,
            exptab$replicate, exptab$direction, MoreArgs = list(allwindchromtib,
            blacklisttib, maptracktib, windsize, currentchrom, chromlength,
            nbcputrans, saveobjectpath, verbose, showtime, showmemory, reload,
            tmpfold), SIMPLIFY = FALSE))
}

.loadbgprocessing <- function(exptab, blacklisttib, maptrackpath, allwindtib,
        windsize, chromtab, nbcputrans, showtime, showmemory, saveobjectpath,
        reload, tmpfold, verbose) {

            if (verbose) message("Removing scores within black list intervals,",
            " keeping those on high mappability regions, and computing ",
            "weighted means.")
            expnamevec <- paste0(exptab$condition, exptab$replicate,
                exptab$direction)

            ## Read the full mappability track once (instead of per chromosome)
            needfullread <- TRUE
            if (reload && !is.na(saveobjectpath)) {
                chromnames <- GenomeInfoDb::seqnames(chromtab)
                cachefiles <- file.path(saveobjectpath,
                    paste0("maptrackbed-", chromnames, ".rds"))
                needfullread <- !all(file.exists(cachefiles))
            }
            fullmaptracktib <- if (needfullread) {
                .readfullmaptrack(maptrackpath, showtime, showmemory, verbose)
            } else NULL

            ## Loading process on a specific chromosom
            invisible(lapply(GenomeInfoDb::seqnames(chromtab),
                function(currentchrom, chromtab, showtime,
                showmemory, saveobjectpath, reload, verbose, exptab,
                blacklisttib, nbcputrans, allwindtib, expnamevec, windsize,
                tmpfold, fullmaptracktib) {

                    if (showtime) start_bglistwmean <- Sys.time()

                    if (verbose) message("\n\t # --- Processing ", currentchrom)
                    ## Filtering the maptrack on the current chromosome
                    chromlength <- .retrievechromlength(chromtab, currentchrom)
                    maptracktib <- .retrievemaptrack(fullmaptracktib, showtime,
                        showmemory, currentchrom, saveobjectpath,
                        reload, verbose)

                    ## Filtering allwindtib on the current chromosome
                    if (verbose) message("\t\t Selecting windows on ",
                        currentchrom)
                    idxchrom <- which(allwindtib$chrom == currentchrom)

                    if (!isTRUE(all.equal(length(idxchrom), 0))) {
                        allwindchromtib <- allwindtib[idxchrom, ]

                        .retrieveandfilterfrombg(exptab, blacklisttib,
                            maptracktib, nbcputrans, allwindchromtib,
                            expnamevec, windsize, currentchrom, chromlength,
                            saveobjectpath, showtime, showmemory, reload,
                            tmpfold, verbose)
                        rm(maptracktib, allwindchromtib)
                        if (showmemory) print(gc()) else invisible(gc())

                        if (showtime) {
                            end_bglistwmean <- Sys.time()
                            timing <- end_bglistwmean - start_bglistwmean
                            message("\t\t ## Built ", currentchrom, " in: ",
                                format(timing, digits = 2))
                        }
                    } else {
                        if (verbose) message("\t\t No transcript annotations ",
                            " found on ", currentchrom, ". Skipping.")
                    }

                }, chromtab, showtime, showmemory, saveobjectpath,
                    reload, verbose, exptab, blacklisttib, nbcputrans,
                    allwindtib, expnamevec, windsize, tmpfold,
                    fullmaptracktib))
}


#' Blacklist High Mappability Regions in Genomic Data
#'
#' @description
#' This function processes genomic data to remove scores that fall within
#' blacklisted regions or have low mappability, and computes weighted means for
#' overlapping windows. The process ensures the integrity of genomic scores by
#' focusing on high mappability regions and excluding blacklisted intervals.
#'
#' @usage
#' blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
#'    nbcputrans, allwindowsbed, windsize, genomename, saveobjectpath = NA,
#'    tmpfold = file.path(tempdir(), "tmptepr"), reload = FALSE, showtime = FALSE,
#'    showmemory = FALSE, chromtab = NA, forcechrom = FALSE, verbose = TRUE)
#'
#' @param maptrackpath Character string. Path to the mappability track file.
#' @param blacklistpath Character string. Path to the blacklist regions file.
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param nbcputrans Number of CPU cores to use for transcript-level operations.
#' @param allwindowsbed Data frame. BED-formatted data frame obtained with the
#'  function 'makewindows'.
#' @param windsize An integer specifying the size of the genomic windows.
#' @param genomename Character string. A valid UCSC genome name. It is used to
#' retrieve chromosome metadata, such as names and lengths.
#' @param saveobjectpath Path to save intermediate R objects. Default is `NA`
#'  and R objects are not saved.
#' @param tmpfold A character string specifying the temporary folder for saving
#'   output files. The temporary files contain the scores for each bedgraph on
#'   each chromosome. Default is \code{file.path(tempdir(), "tmptepr")}.
#' @param reload Logical. If `TRUE`, reloads existing saved objects to avoid
#'  recomputation. Default is `FALSE`. If the function failed during object
#'  saving, make sure to delete the corresponding object.
#' @param showtime A logical value indicating whether to display processing
#'   time.
#' @param showmemory A logical value indicating whether to display memory usage 
#'   during processing.
#' @param chromtab A Seqinfo object containing chromosome information. See
#'  details. Default to NA.
#' @param chromtab A Seqinfo object retrieved with the rtracklayer method
#' SeqinfoForUCSCGenome. If NA, the method is called automatically. Default is
#' NA.
#' @param forcechrom Logical indicating if the presence of non-canonical
#' chromosomes in chromtab (if not NA) should trigger an error. Default is
#' \code{FALSE}.
#' @param verbose A logical value indicating whether to display detailed
#'   processing messages.
#'
#' @return This function does not return a value directly. It saves
#' intermediate results to `tmpfold`. These intermediates files are then
#' combined by the function 'createtablescores'.
#'
#' @details
#' The `blacklisthighmap` function iterates through chromosomes, processes
#' genomic scores by removing those overlapping with blacklisted regions, and
#' ensures that scores within windows are computed using a weighted mean when
#' overlaps occur. The function uses parallel processing for efficiency and
#' supports saving (saveobjectpath) and reloading (reload) intermediate results
#' to optimize workflow.
#'
#' The main steps include:
#' - Reading and processing bedGraph values.
#' - Removing scores overlapping with blacklisted or low mappability regions.
#' - Computing weighted means for overlapping scores in genomic windows.
#' - Saving the processed results to specified path (tmpfold).
#'
#' If chromtab is left to NA, the chromosome information is automatically
#' retrieved from the UCSC server using `genomename`. Otherwise, the Seqinfo
#' object can be retrieved with:
#'      chromtab <- rtracklayer::SeqinfoForUCSCGenome(genomename)
#'
#' @examples
#' \donttest{
#' exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
#' gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
#' maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
#'     package = "tepr")
#' blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
#'     package = "tepr")
#' windsize <- 200
#' genomename <- "hg38"
#' chromtabtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
#' allchromvec <- GenomeInfoDb::seqnames(chromtabtest)
#' chromtabtest <- chromtabtest[allchromvec[which(allchromvec == "chr13")], ]
#'
#' ## Copying bedgraphs to the current directory
#' expdfpre <- read.csv(exptabpath)
#' bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
#'     package = "tepr"))
#' expdfpre$path <- bgpathvec
#' write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
#'     quote = FALSE)
#' exptabpath <- "exptab-preprocessing.csv"
#' 
#' ## Necessary result to call blacklisthighmap
#' allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
#' allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)
#'
#' ## Test blacklisthighmap
#' blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
#'     nbcputrans = 1, allwindowsbed, windsize, genomename,
#'     chromtab = chromtabtest, verbose = FALSE)}
#'
#' @importFrom rtracklayer SeqinfoForUCSCGenome import.bedGraph
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom tibble tibble as_tibble add_column
#' @importFrom dplyr relocate filter
#' @importFrom valr bed_intersect
#' @importFrom methods is
#' @importFrom utils read.csv read.delim write.table
#'
#' @seealso
#' [createtablescores][makewindows]
#'
#' @export


blacklisthighmap <- function(maptrackpath, blacklistpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, genomename = NA, saveobjectpath = NA,
    tmpfold = file.path(tempdir(), "tmptepr"), reload = FALSE, showtime = FALSE,
    showmemory = FALSE, chromtab = NA, forcechrom = FALSE, verbose = TRUE) {

        if (showtime) start_time_fun <- Sys.time()

        if (!file.exists(tmpfold))
            dir.create(tmpfold, recursive = TRUE)

        if (!isTRUE(all.equal(typeof(chromtab), "S4"))) {
            if (is.na(genomename) && is.na(chromtab))
                stop("\n[tepr] Error: Missing genome information.\n",
                    "  Provide 'genomename' or 'chromtab'.\n")

            if (!is.na(chromtab) && !forcechrom) {
                if (!isTRUE(all.equal(is(chromtab), "Seqinfo")))
                    stop("\n[tepr] Error: Invalid chromtab type.\n",
                        "  Must be Seqinfo object from ",
                        "SeqinfoForUCSCGenome().\n")
            }

            ## Retrieving chromosome lengths
            if (is.na(chromtab)) chromtab <- retrievechrom(genomename, verbose)

        } else {
            allchromvec <- GenomeInfoDb::seqnames(chromtab)
            idx <- grep("_|chrM", allchromvec, perl = TRUE, invert = FALSE)
            if (!isTRUE(all.equal(length(idx), 0)))
                stop("\n[tepr] Error: Non-canonical chromosomes in chromtab.\n",
                    "  Set forcechrom=TRUE to proceed anyway.\n")
        }

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- utils::read.csv(exptabpath, header = TRUE)

        if (verbose) message("Reading the black list")
        blacklistbed <- utils::read.delim(blacklistpath, header = FALSE)
        colnames(blacklistbed) <- c("chrom", "start", "end", "type")
        blacklisttib <- tibble::as_tibble(blacklistbed)

        ## Converting allwindowsbed to tib
        colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window")
        allwindtib <- tibble::as_tibble(allwindowsbed)

        ## Cleaning
        rm(blacklistbed, allwindowsbed)
        if (showmemory) print(gc()) else invisible(gc())

        ## Removing scores within black list intervals, keeping those on high
        ## mappability regions, and computing weighted means.
        .loadbgprocessing(exptab, blacklisttib, maptrackpath, allwindtib,
            windsize, chromtab, nbcputrans, showtime, showmemory,
            saveobjectpath, reload, tmpfold, verbose)

        rm(blacklisttib, allwindtib, chromtab)
        if (showmemory) print(gc()) else invisible(gc())

        if (showtime) {
            end_time_fun <- Sys.time()
            timing <- end_time_fun - start_time_fun
            message("\t\t ## All bedgraphs processed in: ",
                format(timing, digits = 2))
        }
}
