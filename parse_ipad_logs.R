library(stringr)
# Note: stringr requires R version 2.11
# Therefore does not currently (15.9.2011) work on CS department computers

veclist <- list() # can i has pointers? no? oh well...

# Finds the previous line with a non-NA value in column col.
# Parameters:
# - curr, the current line in veclist (prev defined in relation to curr)
# - col, column whose previous non-NA value should be found
#
# Returns the index of the previous line with a non-NA value in the given
# column.
find_prev <- function(curr, col) {
    prev = curr - 1
    while (prev > 0 && is.na(veclist[[prev]][[col]]))
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
    while (length(veclist[[nxt]])/2 != max_tps)
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

    for (j in 1:length(lines)) {
        skipi <- c()
        line <- strsplit(lines[j], "//")[[1]][1]
        if (is.na(line)) {
            print("Encountered empty line in the log file.")
            next
        }

        veclist[[j]] <<- lapply(str_extract_all(line, "[0-9]+")[[1]], as.numeric)

        touchpoints <- length(veclist[[j]])/3
        if (touchpoints > max_tps)
            max_tps <- touchpoints;

        # Skip every third value in the vector (these are the "taps: 1" or "taps: 0"
        # values that mostly stay constant) because they do not appear informative.
        for (k in 1:touchpoints)
            skipi <- c(skipi, skip);
        veclist[[j]] <<- veclist[[j]][!skipi]
    }
    return(max_tps);
}

# Fills in empty value pairs in the vector list (global variable veclist)
# with NA values and attempts to parse which features (columns) the values
# belong to when values are missing from the row.
fill_in_empty_values <- function(max_tps) {
    missing_val = c(NA, NA) # what to replace missing (X, Y) value pairs with
    max_seen = FALSE;

    for (j in 1:length(veclist)) {
        if (length(veclist[[j]])/2 == max_tps) # do nothing if no values missing
            next

        prev1 <- find_prev(j, 1)
        prev2 <- find_prev(j, 3)

        if (prev1 == 0) { # missing values are at the beginning of the file
            # need to find the next full row to find out which features the existing
            # values (and missing ones) belong to
            next_max_row <- find_next(max_tps, j);
            
            dist_to_nextx1 <- abs(veclist[[j]][[1]] - veclist[[next_max_row]][[1]])
            dist_to_nextx2 <- abs(veclist[[j]][[1]] - veclist[[next_max_row]][[3]])
            
            # choose the feature with the X value closest to the current X value
            if (dist_to_nextx1 < dist_to_nextx2) {
                veclist[[j]] <<- c(veclist[[j]], missing_val)
            } else
                veclist[[j]] <<- c(missing_val, veclist[[j]])
        }
        else { # missing values are at the middle of end of the file
            dist_to_prevx1 <- abs(veclist[[j]][[1]] - veclist[[prev1]][[1]])
            dist_to_prevx2 <- abs(veclist[[j]][[1]] - veclist[[prev2]][[3]])

            # choose the feature with the X value closest to the current X value
            if (dist_to_prevx1 <= dist_to_prevx2) {
                veclist[[j]] <<- c(veclist[[j]], missing_val)
            }
            else
                veclist[[j]] <<- c(missing_val, veclist[[j]])
        }
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

        result <- do.call(rbind, veclist);
        write.table(result, sprintf("%s/%s", destination, resultfile), sep='\t', row.names=FALSE, col.names=FALSE)
    }
}
