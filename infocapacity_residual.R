# Assume there is a symbolic link called "infocapacity" to infocapacity.R
source("infocapacity")

# Calculates residuals for all sequences in the given directory and writes them
# to files in the result directory (parameter resultdir).
calculate_residuals <- function(directory = ".", resultdir = "residuals", normalize = TRUE)
{
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
        print(nrow(data))
        print(nrow(as.data.frame(residuals)))
        write.table(residuals, sprintf("%s/%s", resultdir, filenames[i]),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE);
    }
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

evaluate_residual_complexity <- function(aligneddir = "aligneddata", residualdir = "residuals", fps = 120)
{
    filenames <- dir(aligneddir, "^[[:digit:]]+_ali_[[:digit:]]+.txt$")
    sequences <- length(filenames)
    all_results <- matrix(nrow = sequences, ncol = 4,
        dimnames = list(1:sequences, c("TP", "RSS", "RSS_resid", "quotient")))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        file1 <- sprintf("%s/%d_ali_%d.txt", aligneddir, j, k)
        file2 <- sprintf("%s/%d_ali_%d.txt", aligneddir, k, j)
        residuals_a <- read.table(sprintf("%s/%02d.txt", residualdir, j))
        residuals_b <- read.table(sprintf("%s/%02d.txt", residualdir, k))

        a <- read.table(file1)
        b <- read.table(file2)

        resid_a <- align_residuals(a, residuals_a)
        resid_b <- align_residuals(b, residuals_b)

        # cut the residual sequences to equal length
        diff <- nrow(resid_a) - nrow(resid_b)
        if (diff < 0)
            resid_b <- resid_b[(abs(diff)+1):nrow(resid_b),]
        if (diff > 0)
            resid_a <- resid_a[(diff+1):nrow(resid_a),]

        data <- remove_duplicate_frames(resid_a, resid_b)
        resid_a <- data[[1]]
        resid_b <- data[[2]]

        n <- nrow(resid_a)
        results_a <- evaluate_residual_shared_information(resid_a, resid_b, n)
        results_b <- evaluate_residual_shared_information(resid_b, resid_a, n)

        all_results <- unpack_results_to_matrix(j, results_a, all_results, fps, nrow(a))
        all_results <- unpack_results_to_matrix(k, results_b, all_results, fps, nrow(b))
    }

    return(all_results)
}
