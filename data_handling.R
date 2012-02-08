# Loads and returns the sequence denoted by the given
# sequence number (assumed to be in the current working directory)
load_sequence <- function(sequence_number) {
    return(read.table(sprintf("%02d.txt", sequence_number)))
}

# Loads the aligned data file "seqnum1_ali_seqnum2.txt" as a data frame
# and returns it.
load_aligned <- function(seqnum1, seqnum2) {
    return(read.table(sprintf("aligneddata/%d_ali_%d.txt", seqnum1, seqnum2)))
}

load_aligned_pair_and_residuals <- function(seqnum1, seqnum2) {
    a <- load_aligned(seqnum1, seqnum2)
    b <- load_aligned(seqnum2, seqnum1)

    residuals_a <- get_aligned_residuals(seqnum1, a)
    residuals_b <- get_aligned_residuals(seqnum2, b)

    residuals <- cut_to_equal_length(residuals_a, residuals_b)
    residuals_a <- residuals[[1]]
    residuals_b <- residuals[[2]]

    return(list(a, residuals_a, b, residuals_b))
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
get_aligned_residuals <- function(sequence_num, aligned_sequence) {
    residuals <- calculate_residuals(sequence_num)
    aligned_residuals <- align_residuals(aligned_sequence, residuals)
    return(aligned_residuals)
}

# Returns the residuals for the given sequence as a data frame.
# Parameters:
# - sequencefile    name of the file containing the coordinate sequence
calculate_residuals <- function(sequence_num) {
    sequence <- load_sequence(sequence_num)
    return(residuals(sequence))
}

# Aligns the given residuals according to the given pre-aligned coordinate
# sequence and returns it in data frame format.
# Parameters:
# - sequence    the aligned coordinate sequence
# - residuals   residuals for the original (non-aligned) sequence
align_residuals <- function(sequence, residuals) {
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
