source("infocapacity.R")

# Returns the residuals for the given sequence as a data frame.
# Parameters:
# - sequencefile    name of the file containing the coordinate sequence
calculate_residuals <- function(sequencefile)
{
    sequence <- normalize_features(read.table(sequencefile))
    n <- nrow(sequence) - 2

    residuals <- evaluate_residuals(cbind(sequence[1:n,],sequence[1:n+1,]), sequence[1:n+2,])

    return(as.data.frame(residuals))
}

# Given a result list, a result matrix and other information used in
# throughput calculation, places the results in the list into their
# places in the matrix and calculates throughput which is also placed
# in the matrix.
# Parameters:
# - sequence_index      index of the sequence, which indicates which row
#                       of the result matrix the results should be placed on
# - results             a result list given by evaluate_residual_shared_information
# - resultmatrix        matrix to unpack the results to
# - fps                 frames per second (to calculate throughut)
# - seq_length          sequence length (after alignment and duplicate removal)
unpack_results_to_matrix <- function(sequence_index, results, resultmatrix, fps,
    seq_length)
{
    quotient <- results$total_RSS / results$total_RSS_residual
    throughput <- results$total_shared / seq_length * fps / log(2.0)

    resultmatrix[sequence_index,1] <- throughput
    resultmatrix[sequence_index,2] <- results$total_RSS
    resultmatrix[sequence_index,3] <- results$total_RSS_residual
    resultmatrix[sequence_index,4] <- quotient
    return(resultmatrix)
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
    residuals <- as.matrix(residuals)
    frames <- nrow(sequence)

    duplicate <- rowSums((sequence[2:frames,]-sequence[1:(frames-1),])^2) == 0
    index <- index_of_third_frame(duplicate)
    duplicate <- duplicate[index:length(duplicate)]

    aligned <- matrix(0, nrow = frames - (index + 1), ncol = ncol(sequence))
    line <- 0

    for (l in 1:(frames - (index + 1))) {
        if (!duplicate[l])
            line <- line + 1

        aligned[l,] <- residuals[line,]
    }

    return(as.data.frame(aligned))
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

# Returns the residuals for a given sequence, aligned according
# to a pre-made residual alignment.
# Parameters:
# - sequencefile     the name of the file containing the original
#                    coordinate sequence
# - aligned_sequence the same sequence after CTW
get_aligned_residuals <- function(sequencefile, aligned_sequence) {
    residuals <- calculate_residuals(sequencefile)
    aligned_residuals <- align_residuals(aligned_sequence, residuals)
    return(aligned_residuals)
}

# Helper function for evaluate_residual complexity.
# Parameters:
# - a, b    aligned residual sequences
pair_residual_complexity <- function(a, b) {
    data <- remove_duplicate_frames(a, b)
    n <- nrow(data[[1]])
    return(evaluate_residual_shared_information(data[[1]], data[[2]], n))
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
# - aligneddir      the name of the subdirectory containing the aligned
#                   coordinate sequences ("aligneddata" by default)
# - fps             frames per second in the given sequences
evaluate_residual_complexity <- function(aligneddir = "aligneddata", fps = 120)
{
    filenames <- dir(aligneddir, "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    sequences <- length(filenames)
    all_results <- matrix(nrow = sequences, ncol = 4,
        dimnames = list(1:sequences, c("TP", "RSS", "RSS_resid", "quotient")))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        a <- read.table(sprintf("%s/%d_ali_%d.txt", aligneddir, j, k))
        b <- read.table(sprintf("%s/%d_ali_%d.txt", aligneddir, k, j))

        residuals_a <- get_aligned_residuals(sprintf("%02d.txt", j), a)
        residuals_b <- get_aligned_residuals(sprintf("%02d.txt", k), b)

        residuals <- cut_to_equal_length(residuals_a, residuals_b)

        results_a <- pair_residual_complexity(residuals[[1]], residuals[[2]])
        results_b <- pair_residual_complexity(residuals[[2]], residuals[[1]])

        n_a <- nrow(residuals[[1]])
        n_b <- nrow(residuals[[2]])

        all_results <- unpack_results_to_matrix(j, results_a, all_results, fps, n_a)
        all_results <- unpack_results_to_matrix(k, results_b, all_results, fps, n_b)
    }

    return(all_results)
}
