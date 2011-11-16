source("infocapacity.R")

calculate_residuals <- function(sequencefile)
{
    sequence <- normalize_features(read.table(sequencefile))
    n <- nrow(sequence) - 2

    residuals <- evaluate_residuals(cbind(sequence[1:n,],sequence[1:n+1,]), sequence[1:n+2,])

    return(as.data.frame(residuals))
}

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

index_of_third_frame <- function(duplicate) {
    startindex <- 0
    for (i in 1:2) {
        startindex <- startindex + 1
        while (duplicate[startindex])
            startindex <- startindex + 1
    }
    return(startindex)
}

align_residuals <- function(sequence, residuals)
{
    residuals <- as.matrix(residuals)
    frames <- nrow(sequence)

    duplicate <- rowSums((sequence[2:frames,]-sequence[1:(frames-1),])^2) < 0.001
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

        residuals_a <- calculate_residuals(sprintf("%02d.txt", j))
        residuals_b <- calculate_residuals(sprintf("%02d.txt", k))

        resid_a <- align_residuals(a, residuals_a)
        resid_b <- align_residuals(b, residuals_b)

        # cut the residual sequences to equal length
        diff <- nrow(resid_a) - nrow(resid_b)
        if (diff < 0)
            resid_b <- resid_b[(abs(diff)+1):nrow(resid_b),]
        if (diff > 0)
            resid_a <- resid_a[(diff+1):nrow(resid_a),]

        data_a <- remove_duplicate_frames(resid_a, resid_b)
        data_b <- remove_duplicate_frames(resid_b, resid_a)

        n_a <- nrow(data_a[[1]])
        n_b <- nrow(data_b[[1]])

        results_a <- evaluate_residual_shared_information(data_a[[1]], data_a[[2]], n_a)
        results_b <- evaluate_residual_shared_information(data_b[[1]], data_b[[2]], n_b)

        all_results <- unpack_results_to_matrix(j, results_a, all_results, fps, nrow(a))
        all_results <- unpack_results_to_matrix(k, results_b, all_results, fps, nrow(b))
    }

    return(all_results)
}
