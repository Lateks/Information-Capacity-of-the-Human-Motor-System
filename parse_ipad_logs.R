library(stringr) # Note: stringr requires R version 2.11

skip = c(FALSE, FALSE, TRUE)
missing_val = c(NA, NA)
max_tps <- 0;
veclist <- list()

find_prev <- function(vecs, curr_i, vec_i) {
    prev = curr_i - 1;
    while (is.na(vecs[[prev]][[vec_i]]))
        prev = prev - 1;
    return(prev);
}

parse_values_to_veclist <- function(lines) {
    for (j in 1:length(lines)) {
        skipi <- c();
        line <- strsplit(lines[j], "//")[[1]][1]
        if (is.na(line))
            next

        veclist[[j]] <<- lapply(str_extract_all(line, "[0-9]+")[[1]], as.numeric)
        touchpoints <- length(veclist[[j]])/3;
        if (touchpoints > max_tps)
            max_tps <<- touchpoints;

        for (k in 1:touchpoints)
            skipi <- c(skipi, skip);
        veclist[[j]] <<- veclist[[j]][!skipi]
    }
}

fill_in_empty_values <- function() {
    max_seen = FALSE;
    for (j in 1:length(veclist)) {
        if (length(veclist[[j]])/2 == max_tps) {
            max_seen <- TRUE
        } else if (!max_seen) {
            veclist[[j]] <<- c(veclist[[j]], missing_val)
        } else { # Note: this works for only 2 touch points
            prev1 <- find_prev(veclist, j, 1)
            prev2 <- find_prev(veclist, j, 3)
            dist_to_prevx1 <- abs(veclist[[j]][[1]] - veclist[[prev1]][[1]])
            dist_to_prevx2 <- abs(veclist[[j]][[1]] - veclist[[prev2]][[3]])
            if (dist_to_prevx1 <= dist_to_prevx2) {
                veclist[[j]] <<- c(veclist[[j]], missing_val)
            }
            else
                veclist[[j]] <<- c(missing_val, veclist[[j]])
        }
    }
}

parse_logs <- function(logdir = ".", destination = "data") {
    logfiles <- dir(logdir, "^[[:digit:]]+.txt$")
    
    for (i in 1:length(logfiles)) {
        logf <- file(sprintf("%s/%s", logdir, logfiles[i]), "r")
        lines <- readLines(logf, -1L)
        close(logf, "r")

        parse_values_to_veclist(lines)

        print(logfiles[i])
        if (max_tps > 2) {
            print("Parsing for more than two touch points not implemented.")
            next
        }

        fill_in_empty_values()

        result <- do.call(rbind, veclist);
        write.table(result, sprintf("%s/%s", destination, logfiles[i]), sep='\t', row.names=FALSE, col.names=FALSE)
    }
}
