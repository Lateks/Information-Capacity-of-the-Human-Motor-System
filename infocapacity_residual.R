# Assume there is a symbolic link called "infocapacity" to infocapacity.R
source("infocapacity")

# Run this first
calculate_residuals <- function(directory = ".", resultdir = "residuals", normalize = TRUE) {
    filenames <- dir(".", "^[[:digit:]]+.txt$")
    if (!file.exists(resultdir))
        dir.create(resultdir)

    for (i in 1:length(filenames)) {
        print(filenames[i])
        data <- read.table(filenames[i])
        if (normalize)
            data <- normalize_features(data)
        n <- nrow(data) - 2

        residuals <- evaluate_residuals(cbind(data[1:n,],data[1:n+1,]), data[1:n+2,])
        write.table(residuals, sprintf("%s/%s", resultdir, filenames[i]),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE);
    }
}

# Then run CTW for the result file pairs.
# Trying out a new, perhaps more sensible(?) naming scheme:
# name all aligned data files x.y_z.txt
# where x    = the number of the sequence "type"
#              (e.g. the glasgow dataset contains 7 sequence types)
#       y, z = repetitions of sequence type x,
#              numbering for each type starts from 1

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

# Run this after you have done the alignment
evaluate_residual_complexity <- function(directory = "residuals/aligneddata", fps = 120)
{
    filenames <- dir(directory, "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    sequences <- length(filenames)
    all_results <- matrix(nrow = sequences, ncol = 4,
        dimnames = list(1:sequences, c("TP", "RSS", "RSS_resid", "quotient")))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        file1 <- sprintf("%s/%d_ali_%d.txt", directory, j, k)
        file2 <- sprintf("%s/%d_ali_%d.txt", directory, k, j)

        a <- read.table(file1)
        b <- read.table(file2)

        data <- remove_duplicate_frames(a, b)
        a <- data[[1]]
        b <- data[[2]]

        results_a <- evaluate_residual_shared_information(a, b)
        results_b <- evaluate_residual_shared_information(b, a)

        all_results <- unpack_results_to_matrix(j, results_a, all_results, fps, nrow(a))
        all_results <- unpack_results_to_matrix(k, results_b, all_results, fps, nrow(b))
    }

    return(all_results)
}

pair_residual_complexity <- function(filename1, filename2, fps = 20)
{
    a <- read.table(filename1)
    b <- read.table(filename2)

    data <- remove_duplicate_frames(a, b)
    a <- data[[1]]
    b <- data[[2]]

    results_a <- evaluate_residual_shared_information(a, b)
    throughput <- results_a$total_shared / nrow(a) * fps / log(2.0)
    return(list(throughput = throughput, results_a))
}
