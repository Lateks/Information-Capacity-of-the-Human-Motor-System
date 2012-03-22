source("infocapacity.R")

# Loads and returns the sequence denoted by the given
# sequence number (assumed to be in the current working directory)
load_sequence <- function(sequence_number) {
    return(read.table(sprintf("%02d.txt", sequence_number)))
}

# Loads the alignment files corresponding to the given
# sequence numbers and calculates residuals.
# Returns a list with the residuals of
# seqnum1_ali_seqnum2 first and the residuals of
# seqnum2_ali_seqnum1 second.
load_aligned_residuals <- function(seqnum1, seqnum2) {
    ali_a <- load_aligned(seqnum1, seqnum2, data = FALSE)
    ali_b <- load_aligned(seqnum2, seqnum1, data = FALSE)

    residuals_a <- get_aligned_residuals(seqnum1, ali_a)
    residuals_b <- get_aligned_residuals(seqnum2, ali_b)

    residuals <- cut_to_equal_length(residuals_a, residuals_b)
    residuals_a <- residuals[[1]]
    residuals_b <- residuals[[2]]

    return(list(residuals_a, residuals_b))
}

# Loads the aligned data file "seqnum1_ali_seqnum2.txt" or the corresponding
# alignment file as a data frame and returns it. (Note: the alignment file
# has only one column, so it is practically a vector.)
#
# Parameters:
# seqnum1       number of the first sequence
# seqnum2       number of the second sequence
# data          logical value: return data (TRUE) or just the alignment (FALSE)
load_aligned <- function(seqnum1, seqnum2, data = TRUE) {
    if (data)
        dir = "aligneddata"
    else
        dir = "alignment"
    return(read.table(sprintf("%s/%d_ali_%d.txt", dir, seqnum1, seqnum2)))
}

# Cuts two sequences to equal length by removing frames from the
# beginning.
#
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
# - alignment        the file containing frame duplication information
#                    in the aligned data (as output by CTW)
get_aligned_residuals <- function(sequence_num, alignment) {
    residuals <- calculate_residuals(sequence_num)
    return(align_residuals(alignment, residuals))
}

# Returns the residuals for the given sequence as a data frame.
# Parameters:
# - sequence_num    the number of the sequence
calculate_residuals <- function(sequence_num) {
    sequence <- load_sequence(sequence_num)
    return(residuals(sequence))
}

# Aligns the given residuals according to the given pre-aligned coordinate
# sequence and returns it in data frame format.
# Parameters:
# - alignment   the alignment file, indicating frame duplications
# - residuals   residuals for the original (non-aligned) sequence
align_residuals <- function(alignment, residuals) {
    frames <- nrow(sequence)
    # The first two frames (numbered 1 and 2) are cut out because
    # there are no residuals for them.
    alignment <- alignment[alignment > 2]  - 2

    aligned <- matrix(0, nrow = length(alignment), ncol = ncol(residuals))
    # Perform alignment according to the given alignment file.
    for (line in 1:nrow(aligned))
        aligned[line,] <- residuals[alignment[line],]

    return(aligned)
}
