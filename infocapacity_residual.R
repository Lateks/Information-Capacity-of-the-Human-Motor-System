source("infocapacity.R")

# Returns the residuals for the given sequence as a data frame.
# Parameters:
# - sequencefile    name of the file containing the coordinate sequence
calculate_residuals <- function(sequence_num)
{
    sequence <- load_sequence(sequence_num)
    n <- nrow(sequence) - 2

    residuals <- evaluate_residuals(cbind(sequence[1:n,],sequence[1:n+1,]), sequence[1:n+2,])

    return(residuals)
}

# Given a result list, a result matrix and other information used in
# throughput calculation, places the results in the list into their
# places in the matrix and calculates throughput which is also placed
# in the matrix.
# Parameters:
# - results             a result list given by evaluate_residual_shared_information
# - fps                 frames per second (to calculate throughut)
# - seq_length          sequence length (after alignment and duplicate removal)
construct_result_vector <- function(results, fps, seq_length)
{
    quotient <- results$RSS / results$RSS_conditional
    throughput <- results$total_shared / seq_length * fps / log(2.0)
    shared_information_bits <- results$total_shared / log(2.0)

    return(c(throughput, shared_information_bits, results$RSS, results$RSS_conditional, quotient))
}

# Returns the index of the third frame of the original sequence
# in the aligned sequence given a logical vector where each
# element indicates whether the corresponding row in the aligned
# sequence is a duplicate of the previous one. The first element
# (indicator for the first row which is never a duplicate) is
# assumed to be missing.
#
# Parameters:
# - duplicate   a logical vector (as described above)
index_of_third_frame <- function(duplicate) {
    startindex <- 0
    for (i in 1:2) {
        startindex <- startindex + 1
        while (duplicate[startindex])
            startindex <- startindex + 1
    }
    return(startindex)
}

# Aligns the given residuals according to the given pre-aligned coordinate
# sequence and returns it in data frame format.
# Parameters:
# - sequence    the aligned coordinate sequence
# - residuals   residuals for the original (non-aligned) sequence
align_residuals <- function(sequence, residuals)
{
    frames <- nrow(sequence)
    features <- ncol(sequence)

    duplicate <- rowSums((sequence[2:frames,]-sequence[1:(frames-1),])^2) == 0
    index <- index_of_third_frame(duplicate)
    duplicate <- duplicate[index:length(duplicate)]

    aligned <- matrix(0, nrow = frames - index, ncol = features)
    line <- 0

    for (l in 1:(frames - index)) {
        if (!duplicate[l])
            line <- line + 1

        aligned[l,] <- residuals[line,]
    }

    return(aligned)
}

# Returns the residuals for a given sequence, aligned according
# to a pre-made residual alignment.
# Parameters:
# - sequencefile     the name of the file containing the original
#                    coordinate sequence
# - aligned_sequence the same sequence after CTW
get_aligned_residuals <- function(sequence_num, aligned_sequence) {
    residuals <- calculate_residuals(sequence_num)
    aligned_residuals <- align_residuals(aligned_sequence, residuals)
    return(aligned_residuals)
}

# Cuts two sequences to equal length by removing frames from the
# beginning.
# Parameters:
# - a, b    residual sequences
cut_to_equal_length <- function(a, b) {
    diff <- nrow(a) - nrow(b)
    if (diff < 0)
        b <- b[(abs(diff)+1):nrow(b),]
    if (diff > 0)
        a <- a[(diff+1):nrow(a),]
    return(list(a, b))
}

# Helper function for evaluate_residual complexity.
# Parameters:
# - a, b    aligned residual sequences
pair_residual_complexity <- function(a, b, fps, pca = FALSE, features = c()) {
    data <- remove_duplicate_frames(a, b)
    a <- data[[1]]
    b <- data[[2]]
    n <- nrow(a)

    if (pca) {
        reduced <- pca(a, b)
        a <- reduced[[1]]
        b <- reduced[[2]]

        for (i in features) {
            if (i != 0) {
                dev.new()
                par(mfcol = c(1, 2))
                plot(a[,i], main = sprintf("Aligned residuals after PCA, feature %d", i))
                plot(b[,i], main = sprintf("Aligned residuals after PCA, feature %d", i))
            }
        }
    }

    results_a <- evaluate_residual_shared_information(a, b)

    if (length(features) > 0 && features[1] == 0) {
        print(results_a$feature_shared / log(2.0))
    }
    else {
        for (i in features)
            print(sprintf("Feature %d shared information in bits: %f", i, results_a$feature_shared[i] / log(2.0)))
    }

    return(construct_result_vector(results_a, fps, n))
}

# Loads and returns the sequence denoted by the given
# sequence number (assumed to be in the current working directory)
load_sequence <- function(sequence_number) {
    return(read.table(sprintf("%02d.txt", sequence_number)))
}

# Plots the features (original and residual, aligned and unaligned)
# given by the vector featurenums.
#
# Parameters:
# - featurenums     vector containing the numbers of the features to plot
# - origa_num       number of the original sequence a
# - origb_num       number of the original sequence b
# - aligneda        sequence a aligned with b
# - alignedb        sequence b aligned with a
# - residuals_a     the aligned residuals for sequence a
# - residuals_b     the aligned residuals for sequence b
plot_features <- function(featurenums, origa_num, origb_num, aligneda, alignedb, residuals_a, residuals_b) {
    if ((length(featurenums) < 1) || (length(featurenums) == 1 && featurenums[1] == 0))
        return()

    origa <- load_sequence(origa_num)
    origb <- load_sequence(origb_num)
    origresiduals_a <- calculate_residuals(origa_num)
    origresiduals_b <- calculate_residuals(origb_num)

    for (i in featurenums) {
        dev.new()
        par(mfcol = c(3, 2))

        # Plot originals
        plot(origa[,i], main = sprintf("Original sequence, feature %d", i))
        points(origb[,i], col = "blue")

        plot(origresiduals_a[,i], main = sprintf("Original residuals, sequence %d, feature %d", origa_num, i))
        plot(origresiduals_b[,i], main = sprintf("Original residuals, sequence %d, feature %d", origb_num, i))

        # Plot aligned
        plot(aligneda[,i], main = sprintf("Aligned sequence, feature %d", i))
        points(alignedb[,i], col = "blue")

        plot(residuals_a[,i], main = sprintf("Aligned residuals, sequence %d, feature %d", origa_num, i))
        plot(residuals_b[,i], main = sprintf("Aligned residuals, sequence %d, feature %d", origb_num, i))
    }
}

# Evaluates residual complexity for a given pair of sequences
# (original sequence files assumed to be in the working directory
# and aligned sequences in the subdirectory "aligneddata").
#
# Parameters:
# - seqnum1, seqnum2    numbers of the sequences
# - fps                 frames per second in the sequences
# - pca                 use PCA (default FALSE)
# - plotfeatures        a vector of feature numbers to plot
#                       (the types of plots given depend on whether
#                       PCA is used or not) and print feature shared
#                       information for (use 0 as the first element
#                       to print feature shared information for all
#                       features instead of just the plotted ones)
evaluate_pair <- function(seqnum1, seqnum2, fps = 120, pca = FALSE, plotfeatures = c())
{
    a <- load_aligned(seqnum1, seqnum2)
    b <- load_aligned(seqnum2, seqnum1)

    residuals_a <- get_aligned_residuals(seqnum1, a)
    residuals_b <- get_aligned_residuals(seqnum2, b)

    residuals <- cut_to_equal_length(residuals_a, residuals_b)
    residuals_a <- residuals[[1]]
    residuals_b <- residuals[[2]]

    if (!pca) {
        plot_features(plotfeatures, seqnum1, seqnum2, a, b, residuals_a, residuals_b)
    }

    results_a <- pair_residual_complexity(residuals_a, residuals_b,
        pca = pca, fps = fps, features = plotfeatures)
    results_b <- pair_residual_complexity(residuals_b, residuals_a,
        pca = pca, fps = fps, features = plotfeatures)

    return(rbind(results_a, results_b, deparse.level = 0))
}

# Calculates residual complexities and throughputs for all sequences in
# the current directory. (Note: the working directory should contain
# the original coordinate sequences. The filenames should be e.g. "01.txt".)
#
# The columns of the result matrix contain throughput, sum of residuals,
# sum of residuals of residuals and the quotient of the previous two,
# in that order.
#
# Parameters:
# - fps             frames per second in the given sequences
# - pca             use principal components analysis
evaluate_residual_complexity <- function(fps = 120, pca = FALSE)
{
    filenames <- dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    sequences <- length(filenames)
    all_results <- matrix(nrow = sequences, ncol = 5,
        dimnames = list(1:sequences, c("TP", "shared", "RSS", "RSS_resid", "quotient")))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        all_results[j:k,] <- evaluate_pair(j, k, fps = fps, pca = pca)
    }

    return(all_results)
}

# Calculate the number of k-combinations in a set with n elements.
combinations <- function(n, k) {
    return(factorial(n) / (factorial(k) * factorial(n - k)))
}

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
# - compare boolean value indicating that the throughput results
#           from earlier methods should be included
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
subdir_based_residual_complexity <- function(fps = 120, pca = FALSE, compare = FALSE) {
    subdirs <- dir(".", "^[[:digit:]][[:digit:]]+$")
    results <- list()
    for (i in 1:length(subdirs)) {
        setwd(subdirs[i])
        sequences <- length(dir(".", "^[[:digit:]]+.txt$"))

        rownames <- c()
        combos <- combinations(sequences, 2) * 2
        col_names <- c("TP", "shared", "RSS", "RSS_resid", "quotient")
        all_results <- matrix(nrow = combos, ncol = 5,
            dimnames = list(1:combos, col_names))
        rows <- c(1, 2)

        if (compare)
            compare_results <- matrix(nrow = combos, ncol = 2,
                dimnames = list(1:combos, c("original_TP", "old_residual_TP")))

        for (j in 1:(sequences-1)) {
            for (k in (j+1):sequences) {
                rownames <- c(rownames, sprintf("(%d, %d)", j, k))
                rownames <- c(rownames, sprintf("(%d, %d)", k, j))

                all_results[rows[1]:rows[2],] <- evaluate_pair(j, k, fps = fps, pca = pca)

                if (compare) { # run earlier method versions for the same files
                    result1 <- pair_throughput(j, k, fps = fps, pca = pca, residuals = FALSE)
                    result2 <- pair_throughput(k, j, fps = fps, pca = pca, residuals = FALSE)
                    compare_results[rows[1],] <- c(result1$throughput, result2$throughput)

                    result3 <- pair_throughput(j, k, fps = fps, pca = pca)
                    result4 <- pair_throughput(k, j, fps = fps, pca = pca)
                    compare_results[rows[2],] <- c(result3$throughput, result4$throughput)
                }
                rows <- rows + 2
            }
        }
        rownames(all_results) <- rownames
        if (compare)
            all_results <- cbind(all_results, compare_results)

        results[[i]] <- list(subdirs[i], all_results);
        setwd("..")
    }
    return(results)
}
