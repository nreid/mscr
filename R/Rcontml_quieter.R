#this function (and the relevant helper function) are taken from package Rphylip by Liam J. Revell and Scott A. Chamberlain
#I modified it to be quieter. It requires a helper function that isn't exported, and so can't be seen by this function, so I included it. 

Rcontml_quieter <- function (X, path = NULL, ...) 
{
    if (is.null(path)) 
        path <- findPath("contml")
    if (is.null(path)) 
        stop("No path provided and was not able to find path to contml")
    if (hasArg(quiet)) 
        quiet <- list(...)$quiet
    else quiet <- FALSE
    if (!quiet) 
        if (file.warn(c("infile", "outfile", "outtree")) == 0) 
            return(NULL)
    oo <- c("r")
    if (is.matrix(X)) {
        write(paste("    ", nrow(X), "   ", ncol(X), sep = ""), 
            file = "infile")
        for (i in 1:nrow(X)) {
            sp <- as.character(i)
            sp <- paste(sp, paste(rep(" ", 11 - nchar(sp)), collapse = ""), 
                collapse = "")
            tt <- paste(sp, paste(X[i, ], collapse = " "), collapse = " ")
            write(tt, append = TRUE, file = "infile")
        }
        oo <- c(oo, "c")
        if (hasArg(tree)) {
            oo <- c(oo, "u")
            tree <- list(...)$tree
            tree$tip.label <- sapply(tree$tip.label, function(x, 
                y) which(x == y), y = rownames(X))
            write.tree(tree, "intree")
            intree <- TRUE
        }
        else intree <- FALSE
        if (hasArg(global)) 
            global <- list(...)$global
        else global <- TRUE
        if (global) 
            oo <- c(oo, "g")
        if (hasArg(random.order)) 
            random.order <- list(...)$random.order
        else random.order <- TRUE
        if (random.order) {
            if (hasArg(random.addition)) 
                random.addition <- list(...)$random.addition
            else random.addition <- 10
            oo <- c(oo, "j", sample(seq(1, 99999, by = 2), 1), 
                random.addition)
        }
        if (quiet) 
            oo <- c(oo, 2)
        oo <- c(oo, "y", "r")
        system("touch outfile")
        system(paste(path, "/contml", sep = ""), input = oo, 
            ignore.stdout = TRUE, ignore.stderr = TRUE)
        tree <- read.tree("outtree")
        temp <- readLines("outfile")
        logLik <- as.numeric(strsplit(temp[grep("Ln Likelihood", 
            temp)], "=")[[1]][2])
        if (!quiet){
            temp <- lapply(temp, function(x) {
                cat(x)
                cat("\n")
                })
            }          
        if (!quiet) {
            cat("Translation table\n")
            cat("-----------------\n")
            temp <- lapply(1:nrow(X), function(x, y) cat(paste("\t", 
                paste(x, y[x], sep = "\t"), "\n", sep = "")), 
                y = rownames(X))
            cat("\n")
        }
        tree$tip.label <- rownames(X)[as.numeric(tree$tip.label)]
    }
    else if (is.list(X)) {
        tips <- rownames(X[[1]])
        X <- lapply(X, function(x, tips) x[tips, ], tips = tips)
        write(paste("    ", nrow(X[[1]]), "   ", length(X), sep = ""), 
            file = "infile")
        nalleles <- sapply(X, ncol)
        write(paste(nalleles, collapse = " "), file = "infile", 
            append = TRUE)
        temp <- sapply(X, rowSums)
        if (!all(round(temp, 2) == 1)) 
            stop("Some of the rows of X do not sum to 1.0")
        for (i in 1:length(tips)) {
            sp <- as.character(i)
            sp <- paste(sp, paste(rep(" ", 11 - nchar(sp)), collapse = ""), 
                collapse = "")
            dd <- vector()
            for (j in 1:length(X)) dd <- c(dd, X[[j]][i, ])
            tt <- paste(sp, paste(dd, collapse = " "), collapse = " ")
            write(tt, append = TRUE, file = "infile")
        }
        oo <- c(oo, "a")
        if (hasArg(tree)) {
            oo <- c(oo, "u")
            tree <- list(...)$tree
            tree$tip.label <- sapply(tree$tip.label, function(x, 
                y) which(x == y), y = rownames(X))
            write.tree(tree, "intree")
            intree <- TRUE
        }
        else intree <- FALSE
        if (hasArg(global)) 
            global <- list(...)$global
        else global <- TRUE
        if (global) 
            oo <- c(oo, "g")
        if (hasArg(random.order)) 
            random.order <- list(...)$random.order
        else random.order <- TRUE
        if (random.order) {
            if (hasArg(random.addition)) 
                random.addition <- list(...)$random.addition
            else random.addition <- 10
            oo <- c(oo, "j", sample(seq(1, 99999, by = 2), 1), 
                random.addition)
        }
        if (quiet) 
            oo <- c(oo, 2)
        oo <- c(oo, "y", "r")
        system("touch outfile")
        system(paste(path, "/contml", sep = ""), input = oo, 
            ignore.stdout = TRUE, ignore.stderr = TRUE)
        tree <- read.tree("outtree")
        temp <- readLines("outfile")
        logLik <- as.numeric(strsplit(temp[grep("Ln Likelihood", 
            temp)], "=")[[1]][2])
        if (!quiet) {
            temp <- lapply(temp, function(x) {
                cat(x)
                cat("\n")
                })
            }
        if (!quiet) {
            cat("Translation table\n")
            cat("-----------------\n")
            temp <- lapply(1:length(tips), function(x, y) cat(paste("\t", 
                paste(x, y[x], sep = "\t"), "\n", sep = "")), 
                y = tips)
            cat("\n")
        }
        tree$tip.label <- tips[as.numeric(tree$tip.label)]
    }
    else stop("X should be a matrix (for continuous characters) or a list (for gene frequencies)")
    if (hasArg(outgroup)) {
        outgroup <- list(...)$outgroup
        tree <- outgroup.root(tree, outgroup, quiet)
    }
    if (hasArg(cleanup)) 
        cleanup <- list(...)$cleanup
    else cleanup <- TRUE
    if (cleanup) {
        files <- c("infile", "outfile", "outtree")
        if (intree) 
            files <- c(files, "intree")
        cleanFiles(files)
    }
    tree$logLik <- logLik
    return(tree)
}


cleanFiles<-function (fs) 
{
    if (.Platform$OS.type == "windows") 
        for (i in 1:length(fs)) system(paste("rm", fs[i], sep = " "), 
            show.output.on.console = FALSE)
    else for (i in 1:length(fs)) system(paste("rm", fs[i], sep = " "))
}
