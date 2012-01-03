source("infocapacity_residual.R")

# Returns a logical matrix of the same size as the input
# matrix residuals. Each element is either TRUE (the
# corresponding element in the residuals matrix is
# considered an outlier) or FALSE (it is not an outlier).
#
# An "outlier" is defined to be a residual observation that
# differs more than limit (a parameter) times the standard
# deviation of the feature from the feature mean.
#
# Parameters:
# - residuals       residual matrix
# - limit           see definition of an outlier above
detect_outliers <- function(residuals, limit) {
    num_features <- ncol(residuals)
    outliers <- matrix(data = FALSE, nrow = nrow(residuals), ncol = num_features)

    for (i in 1:num_features) {
        sdev <- sd(residuals[,i])
        feat_mean <- mean(residuals[,i])
        outliers[,i] <- abs(feat_mean - residuals[,i]) > limit * sdev
    }

    return(outliers)
}

# Calculates the number of "outliers" in each feature from
# a sequence of residuals and throughput for each of the
# features. Results are returned in a features x 2 matrix
# where the first column contains the number of outliers
# in each feature and the second column contains throughput
# for the same feature.
#
# Parameters:
# - main_sequence_number        number of the sequence to calculate
#                               outliers for
# - sideinfo_sequence_number    number of the sequence to use as
#                               side information for throughput
#                               calculation
# - limit                       see definition of an outlier in
#                               detect_outliers
# - fps                         frames per second in the sequences
#                               (for throughput calculation)
count_outliers <- function(main_sequence_number, sideinfo_sequence_number, limit, fps) {
    residuals <- calculate_residuals(main_sequence_number)
    num_features <- ncol(residuals)
    outliers <- vector(mode = "integer", length = num_features)
    detected <- detect_outliers(residuals, limit)

    for (i in 1:num_features) {
        sdev <- sd(residuals[,i])
        feat_mean <- mean(residuals[,i])
        outliers[i] <- sum(detected[,i])
    }

    data <- load_aligned_pair_and_residuals(main_sequence_number,
        sideinfo_sequence_number)
    residual_seqs <- remove_duplicate_frames(data[[2]], data[[4]])
    feature_tps <- (evaluate_residual_shared_information(residual_seqs[[1]],
        residual_seqs[[2]]))$feature_shared / nrow(residual_seqs[[1]]) * fps / log(2.0);

    outliers_and_throughputs <- cbind(outliers, feature_tps)

    return(outliers_and_throughputs)
}

# Tests different limits in the range specified by the limit_range
# parameter. (Default from 5 to 20.) Interval size is 1.
test_limits <- function(limit_range = c(5, 20), fps = 120) {
    files <- length(dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$"))
    par(mfrow = c(ceiling((limit_range[2] - limit_range[1] + 1)/5), 5))

    for (limit in seq(limit_range[1], limit_range[2], 1)) {
        all_outliers <- c()
        for (i in 1:(files/2)) {
            j = 2 * i - 1
            k = j + 1

            outliers <- count_outliers(j, k, limit, fps)

            all_outliers <- rbind(all_outliers, outliers)
        }
        plot_model(all_outliers, limit)
    }
}

# Plots the outlier data, fits a linear model into it
# (also plotted) and runs ANOVA on the model (results
# printed into the R terminal).
plot_model <- function(all_outliers, limit) {
    max_outliers <- max(all_outliers[,1])
    max_tp <- floor(max(all_outliers[,2]))
    min_tp <- floor(min(all_outliers[,2]))

    outlierfit <- aov(all_outliers[,2] ~ all_outliers[,1])

    # Plot number of outliers vs throughput and a
    # linear model that has been fitted to these values.
    plot(all_outliers, xlim = c(0, max_outliers), ylim = c(min_tp, max_tp),
        main = "Outliers vs. feature throughput",
        sub = sprintf("Limit = %d", limit),
        xlab = "Number of outliers", ylab = "Throughput")
    abline(outlierfit, col = "red")

    print(limit)
    print(anova(outlierfit))
}

# Given a matrix containing outlier counts (first column) and
# feature throughputs (second column), return a matrix with
# feature throughput means (first column) and standard deviations
# (second column) for a particular outlier count (indicated by
# the row number in the result matrix). The number of occurrences
# of a certain outlier count is returned in the third column.
means_and_sdevs <- function(all_outliers) {
    max_outliers <- max(all_outliers[,1])
    tp_means <- matrix(data = 0, nrow = max_outliers+1, ncol = 3,
        dimnames = list(c(0:max_outliers), c("mean tp", "sdev", "occurrences")))
    for (i in 0:max_outliers) {
        ivals = all_outliers[all_outliers[,1] == i,2]
        tp_means[i+1, 1] <- mean(ivals)
        tp_means[i+1, 2] <- sd(ivals)
        tp_means[i+1, 3] <- length(ivals)
    }
    return(tp_means)
}

# Draws useful plots of the outliers for sequences
# in the aligneddata subdirectory of the working
# directory.
#
# Parameters:
# - limit               limit for detecting an outlier
#                       (see detect_outliers)
# - fps                 framerate
plot_outliers <- function(limit = 7, fps = 120) {
    files <- length(dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$"))
    all_outliers <- c()

    for (i in 1:(files/2)) {
        j = 2 * i - 1
        k = j + 1

        outliers <- count_outliers(j, k, limit, fps)

        all_outliers <- rbind(all_outliers, outliers)
    }

    max_outliers <- max(all_outliers[,1])
    max_tp <- floor(max(all_outliers[,2]))
    min_tp <- floor(min(all_outliers[,2]))

    par(mfrow = c(2, 2))
    plot_model(all_outliers, limit)

    tp_means <- means_and_sdevs(all_outliers)

    # Plot the number of occurrences of each number
    # of outliers.
    barplot(tp_means[,3], names.arg = 0:max_outliers,
        main = "Number of occurrences",
        ylab = "Occurrences", xlab = "Number of outliers")

    # Plot throughput means and standard deviations against
    # the number of outliers.
    plot(0:max_outliers, tp_means[,1], col = "red", pch = 18,
        main = "Throughput mean", xlab = "Number of outliers",
        ylab = "Throughput mean")
    plot(0:max_outliers, tp_means[,2], col = "blue", pch = 18,
        main = "Throughput standard deviation",
        xlab = "Number of outliers", ylab = "Throughput sdev")
}

# Plot aligned residuals for the two given sequences (a pair). Only
# the features specified in the features list are plotted (an
# individual window will open for each feature).
#
# Parameters:
# - sequence_number     the number of the first sequence
# - pair_number         the number of the pair of the first sequence
# - features            a vector of the feature numbers to be plotted
# - limit               limit for detecting an outlier
#                       (see detect_outliers)
# - fps                 framerate
show_outliers <- function(sequence_number, pair_number, features = c(1),
    limit = 7, fps = 120) {

    aligned <- load_aligned_pair_and_residuals(sequence_number, pair_number)
    a <- aligned[[1]]
    b <- aligned[[3]]
    residuals_a <- aligned[[2]]
    residuals_b <- aligned[[4]]

    frames <- nrow(residuals_a)
    framenums <- c(1:frames)

    outliers_a <- detect_outliers(residuals_a, limit)
    outliers_b <- detect_outliers(residuals_b, limit)

    for (feat in features) {
        dev.new()
        par(mfrow = c(1, 2))

        feat_outliers_a <- residuals_a[,feat][outliers_a[,feat]]
        if (length(feat_outliers_a) > 0) {
            plot(framenums[outliers_a[,feat]], feat_outliers_a, col = "red",
                xlim = c(1, frames), main = sprintf("Sequence %d, feature %d",
                sequence_number, feat), xlab = "Frames", ylab = "Residuals",
                ylim = c(min(residuals_a[,feat]), max(residuals_a[,feat])))
            points(framenums[!outliers_a[,feat]], residuals_a[,feat][!outliers_a[,feat]])
        }
        else {
            plot(residuals_a[,feat], xlim = c(1, frames),
                main = sprintf("Sequence %d, feature %d", sequence_number,
                feat), xlab = "Frames", ylab = "Residuals")
        }

        feat_outliers_b <- residuals_b[,feat][outliers_b[,feat]]
        if (length(feat_outliers_b) > 0) {
            plot(framenums[outliers_b[,feat]], feat_outliers_b, col = "red",
                xlim = c(1, frames), main = sprintf("Sequence %d, feature %d",
                pair_number, feat), xlab = "Frames", ylab = "Residuals",
                ylim = c(min(residuals_b[,feat]), max(residuals_b[,feat])))
            points(framenums[!outliers_b[,feat]], residuals_b[,feat][!outliers_b[,feat]])
        }
        else {
            plot(residuals_b[,feat], xlim = c(1, frames),
                main = sprintf("Sequence %d, feature %d", pair_number, feat),
                xlab = "Frames", ylab = "Residuals")
        }
    }
}
