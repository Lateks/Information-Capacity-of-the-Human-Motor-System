# Parse iPad multitouch logs into matrices usable by the infocapacity
# script. Requires R version >= 2.11.
#
# Load the file in R and call the parse_logs(). The result files are
# saved into a subdirectory which is by default called "data". The
# current directory is assumed to contain the log files (named
# <number>.log).
#
# The parse_logs function also optionally takes two parameters: logdir
# which is the directory containing the logs and destination which is the
# directory the results will be written to.
#
# Note that parsing ipad logs reliably is practically impossible in
# some cases due to the lack of identifiers on touch points. Therefore
# the results should always be checked and if necessary, corrected by
# hand. (You can also change the guess_correct_features function slightly
# to look at Y coordinates instead of X coordinates if they seem like
# better indicators in your data. Brief instructions included in the function
# comments.)
#
# Also note that parsing for more than two simultaneous touch points is
# not implemented in the current version due to lack of test data. (Feel
# free to add this functionality if required.)
# 
# ======
# TL;DR:
# Usage: parse_logs() or
# parse_logs("ipad_logfile_directory", "result_directory")
# Check results before further use.
# Does not work with data that has more than two simultaneous touch
# points. Support for more touch points should be implemented.
# ======
#
# by Laura Leppänen, 21st September 2011

library(stringr)
# Note: stringr requires R version 2.11
# Therefore does not currently (15.9.2011) work on CS department computers

veclist <- list() # can i has pointers? no? oh well...

# Finds the previous line with a non-NaN value in column col.
# Parameters:
# - curr, the current line in veclist (prev defined in relation to curr)
# - col, column whose previous non-NaN value should be found
#
# Returns the index of the previous line with a non-NaN value in the given
# column or zero if none are found.
find_prev <- function(curr, col) {
    prev = curr - 1
    while (prev > 0 && is.nan(veclist[[prev]][[col]]))
        prev = prev - 1
    return(prev)
}

# Finds the next full line (vector) in veclist. Parameters:
# - max_tps, maximum number of simultaneous touch points (needed to identify
#            full lines)
# - curr, the current line in veclist (next defined in relation to curr)
# 
# Returns the index of the next full line.
find_next <- function(max_tps, curr) {
    nxt = curr + 1
    while ((length(veclist[[nxt]]) - 1)/2 != max_tps)
        nxt = nxt + 1
    return(nxt)
}

# Parses integer values from the given list of lines (a character vector)
# into the global variable veclist. Vectors can be of different lengths
# if some values are missing on some lines.
#
# Returns the maximum number of simultaneous touch points in the log file.
parse_values_to_veclist <- function(lines) {
    max_tps <- 0
    skip = c(FALSE, FALSE, TRUE)
    last_timestamp = list(0.0)

    for (j in 1:length(lines)) {
        skipi <- c(FALSE) # the first value (ms) will never be skipped
        line <- strsplit(lines[j], "//")[[1]]
        if (is.na(line[1])) {
            print("Encountered empty line in the log file.")
            next
        }

        time <- list(as.numeric(str_extract(line[2], "[0-9]+\\.[0-9]+")))
        pointdata <- lapply(str_extract_all(line[1], "[0-9]+")[[1]], as.numeric)

        timestamp <- list(last_timestamp[[1]] + time[[1]])
        veclist[[j]] <<- c(timestamp, pointdata)
        last_timestamp <- timestamp

        touchpoints <- (length(veclist[[j]]) - 1)/3
        if (touchpoints > max_tps)
            max_tps <- touchpoints

        # Skip every third value in the vector/list (these are the "taps: 1" or
        # "taps: 0" values that mostly stay constant) because they do not appear
        # informative.
        for (k in 1:touchpoints)
            skipi <- c(skipi, skip)
        veclist[[j]] <<- veclist[[j]][!skipi]
    }
    return(max_tps)
}

# Compares the distances of the X variable on the current row (curr)
# in veclist to the X variables on rows other1 (first X variable) and
# other2 (second X variable).
#
# Assume that other2=other1 when other2 is not explicitly given.
#
# Estimates which feature the X variable values on the current row
# belong to according to distance comparison results.
#
# Finally returns the current row vector with missing values replaced
# with NaN.
guess_correct_features <- function(curr, other1, other2 = NULL) {
    missing_val <- c(NaN, NaN) # what to replace missing (X, Y) value pairs with
    if (is.null(other2))
        other2 <- other1

    dist_to_otherx1 <- abs(veclist[[curr]][[2]] - veclist[[other1]][[2]])
    dist_to_otherx2 <- abs(veclist[[curr]][[2]] - veclist[[other2]][[4]])

    # Y coordinates could be used instead of X coordinates if they produce
    # better results in some circumstances. Remember to change the if
    # statement below.
    #dist_to_othery1 <- abs(veclist[[curr]][[3]] - veclist[[other1]][[3]])
    #dist_to_othery2 <- abs(veclist[[curr]][[3]] - veclist[[other2]][[5]])
    
    # Choose the X feature with values closest to the values on the current
    # row. A word of caution: this might not work so wellwith very fast
    # touch point movements.
    if (dist_to_otherx1 < dist_to_otherx2) {
        return(c(veclist[[curr]], missing_val))
    } else
        return(c(veclist[[curr]][1], missing_val, veclist[[curr]][2:3]))
}

# Fills in empty value pairs in the vector list (global variable veclist)
# with NaN values and attempts to parse which features (columns) the values
# belong to when values are missing from the row.
fill_in_empty_values <- function(max_tps) {
    max_seen = FALSE

    for (j in 1:length(veclist)) {
        if ((length(veclist[[j]]) - 1)/2 == max_tps) # do nothing if no values missing
            next

        prev1 <- find_prev(j, 2)
        prev2 <- find_prev(j, 4)

        if (prev1 == 0) { # missing values are at the beginning of the file
            # need to find the next full row to find out which features the existing
            # values (and missing ones) belong to
            next_max_row <- find_next(max_tps, j)
            veclist[[j]] <<- guess_correct_features(j, next_max_row)
        } else # missing values are at the middle or end of the file
            veclist[[j]] <<- guess_correct_features(j, prev1, prev2)
    }
}

# This is the main function that should be called. Parameters are:
# - logdir, path to the directory containing the multitouch logs ("." by default)
# - destination, path to the directory that should contain the parsed data ("data" by default)
#
# It is assumed that log files are named <number>.log (e.g. "15.log"). Corresponding result
# data files will be named <number>.txt.
parse_logs <- function(logdir = ".", destination = "data") {
    logfiles <- dir(logdir, "^[[:digit:]]+.log$")
    if (!file.exists(destination))
        dir.create(destination)
    
    for (i in 1:length(logfiles)) {
        print(logfiles[i])
        veclist <<- list()
        resultfile = sprintf("%s.txt", strsplit(logfiles[i], "[[:punct:]]")[[1]][1])

        logf <- file(sprintf("%s/%s", logdir, logfiles[i]), "r")
        lines <- readLines(logf, -1L)
        close(logf, "r")

        max_tps <- parse_values_to_veclist(lines)

        if (max_tps > 2) {
            print("Parsing for more than two touch points not implemented.")
            next
        }

        fill_in_empty_values(max_tps)

        result <- do.call(rbind, veclist)
        write.table(result, sprintf("%s/%s", destination, resultfile), sep='\t', row.names=FALSE, col.names=FALSE)
    }
}
