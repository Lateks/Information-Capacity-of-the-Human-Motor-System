source("infocapacity_residual.R")

# Calculates the number of "outliers" in each feature from
# a sequence of residuals and throughput for each of the
# features. Results are returned in a features x 2 matrix
# where the first column contains the number of outliers
# in each feature and the second column contains throughput
# for the same feature.
#
# An "outlier" is defined to be a residual observation that
# differs more than limit (a parameter) times the standard
# deviation of the feature from the feature mean.
#
# Parameters:
# - main_sequence_number        number of the sequence to calculate
#                               outliers for
# - sideinfo_sequence_number    number of the sequence to use as
#                               side information for throughput
#                               calculation
# - limit                       see definition of an outlier above
# - fps                         frames per second in the sequences
#                               (for throughput calculation)
detect_outliers <- function(main_sequence_number, sideinfo_sequence_number, limit, fps) {
    residuals <- calculate_residuals(main_sequence_number)
    num_features <- ncol(residuals)
    outliers <- vector(mode = "integer", length = num_features)

    for (i in 1:num_features) {
        sdev <- sd(residuals[,i])
        feat_mean <- mean(residuals[,i])
        outliers[i] <- length(residuals[abs(feat_mean - residuals[,i]) > limit * sdev, i])
    }

    data <- load_aligned_pair_and_residuals(main_sequence_number, sideinfo_sequence_number)
    residual_seqs <- remove_duplicate_frames(data[[2]], data[[4]])
    feature_tps <- (evaluate_residual_shared_information(residual_seqs[[1]], residual_seqs[[2]]))$feature_shared / nrow(residual_seqs[[1]]) * fps / log(2.0);

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

            outliers1 <- detect_outliers(j, k, limit, fps)

            all_outliers <- rbind(all_outliers, outliers1)
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
plot_outliers <- function(limit = 6, fps = 120) {
    files <- length(dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$"))
    all_outliers <- c()

    for (i in 1:(files/2)) {
        j = 2 * i - 1
        k = j + 1

        outliers <- detect_outliers(j, k, limit, fps)

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