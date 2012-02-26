source("infocapacity.R")
source("data_handling.R")
library("methods")

residual_complexity <- function(fps = 120, pca = FALSE) {
    filenames <- dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    sequences <- length(filenames)
    all_results <- matrix(nrow = sequences, ncol = 5,
        dimnames = list(1:sequences, c("TP", "shared", "RSS", "RSS_resid", "quotient")))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        all_results[j:k,] <- pair_residual_complexity(j, k, fps = fps, pca = pca)
    }

    return(all_results)
}

pair_residual_complexity <- function(seqnum_a, seqnum_b, fps = 120, pca = FALSE) {
    pair_a_b <- new("residualComplexityEvaluator", seqnum_a = seqnum_a,
                    seqnum_b = seqnum_b, fps = fps, pca = pca)
    results_a_b <- evaluate_complexity(pair_a_b)
    pair_b_a <- new("residualComplexityEvaluator", seqnum_a = seqnum_b,
                    seqnum_b = seqnum_a, fps = fps, pca = pca)
    results_b_a <- evaluate_complexity(pair_b_a)

    return(rbind(results_a_b, results_b_a, deparse.level = 0))
}

setClass("residualComplexityEvaluator",
         representation(seqnum_a = "numeric",
                        seqnum_b = "numeric",
                        fps   = "numeric",
                        pca   = "logical",
                        seq_a = "matrix",
                        seq_b = "matrix"))

setGeneric("evaluate_complexity",
           function(this) standardGeneric("evaluate_complexity"))
setMethod("evaluate_complexity", "residualComplexityEvaluator",
           function(this) {
               this <- load_residuals(this)
               this <- remove_duplicates(this)
               if (this@pca) this <- do_pca(this)
               results <- evaluate_residual_shared_information(this@seq_a, this@seq_b)
               construct_result_vector(this, results)
           })

setGeneric("load_residuals",
           function (this) standardGeneric("load_residuals"))
setMethod("load_residuals", "residualComplexityEvaluator",
          function(this) {
              data <- load_aligned_residuals(this@seqnum_a, this@seqnum_b)
              this@seq_a <- data[[1]]
              this@seq_b <- data[[2]]
              this
          })

setGeneric("remove_duplicates",
           function(this) standardGeneric("remove_duplicates"))
setMethod("remove_duplicates", "residualComplexityEvaluator",
          function(this) {
              data <- remove_duplicate_frames(this@seq_a, this@seq_b)
              this@seq_a <- data[[1]]
              this@seq_b <- data[[2]]
              this
          })

setGeneric("do_pca",
           function(this) standardGeneric("do_pca"))
setMethod("do_pca", "residualComplexityEvaluator",
          function(this) {
              eigenvectors <- pca(this@seq_a)
              this@seq_a <- normalize_features(this@seq_a) %*% eigenvectors
              this@seq_b <- normalize_features(this@seq_b) %*% eigenvectors
              this
          })

# Generates a vector with throughput, shared information
# in bits, total sum of residuals, total conditional sum of
# residuals and the residual quotient. (Vector length is 5.)
setGeneric("construct_result_vector",
           function(this, results) standardGeneric("construct_result_vector"))
setMethod("construct_result_vector", "residualComplexityEvaluator",
          function(this, results) {
              quotient <- results$RSS / results$RSS_conditional
              throughput <- results$total_shared / nrow(this@seq_a) * this@fps / log(2.0)
              shared_information_bits <- results$total_shared / log(2.0)

              c(throughput, shared_information_bits, results$RSS,
                results$RSS_conditional, quotient)
          })

# Calculate residual complexity assuming that each sequence type is in
# its own subdirectory. The original sequence files for different
# sequences should be in subdirectories named 01, 02, 03 etc.
# These subdirectories in turn should all have a subdirectory called
# aligneddata that contains all sequence repetitions aligned with
# each other. See below for an example.
#
# Parameters:
# - fps     frames per second in the sequences
# - pca     use PCA (default FALSE)
#
# Example of the directory tree:
# - working directory
#   |
#   -- 01 -- contains three repetitions of sequence #1
#   |  |
#   |  * 01.txt
#   |  * 02.txt
#   |  * 03.txt
#   |  -- aligneddata
#   |     |
#   |     * 1_ali_2.txt
#   |     * 1_ali_3.txt
#   |     * ...
#   -- 02 -- contains two repetitions of sequence #2
#      |
#      * 01.txt
#      * 02.txt
#      -- aligneddata
#         |
#         * 1_ali_2.txt
#         * 2_ali_1.txt
subdir_based_residual_complexity <- function(fps = 120, pca = FALSE) {
    subdirs <- dir(".", "^[[:digit:]][[:digit:]]+$")
    results <- list()
    for (i in 1:length(subdirs)) {
        setwd(subdirs[i])
        sequences <- length(dir(".", "^[[:digit:]]+.txt$"))

        rownames <- c()
        combos <- combinations(sequences, 2) * 2
        col_names <- c("TP", "shared", "RSS", "RSS_cond", "quotient")
        all_results <- matrix(nrow = combos, ncol = 5,
            dimnames = list(1:combos, col_names))
        rows <- c(1, 2)

        for (j in 1:(sequences-1)) {
            for (k in (j+1):sequences) {
                rownames = append(rownames, sprintf("(%d, %d)", j, k))
                rownames = append(rownames, sprintf("(%d, %d)", k, j))

                all_results[rows[1]:rows[2],] <- pair_residual_complexity(j, k, fps, pca)

                rows <- rows + 2
            }
        }
        rownames(all_results) <- rownames
        results[[i]] <- list(dir = subdirs[i], results = all_results);

        setwd("..")
    }
    return(results)
}

# Calculate the number of k-combinations in a set with n elements.
combinations <- function(n, k) {
    return(factorial(n) / (factorial(k) * factorial(n - k)))
}
