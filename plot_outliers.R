source("infocapacity_residual.R")

detect_outliers <- function(main_sequence_number, sideinfo_sequence_number, limit, fps) {
    residuals <- calculate_residuals(main_sequence_number)
    num_features <- ncol(residuals)
    outliers <- vector(mode = "integer", length = num_features)

    for (i in 1:num_features) {
        sdev <- sd(residuals[,i])
        outliers[i] <- length(residuals[abs(residuals[,i]) > limit * sdev, i])
    }

    data <- load_aligned_pair_and_residuals(main_sequence_number, sideinfo_sequence_number)
    residual_seqs <- remove_duplicate_frames(data[[2]], data[[4]])
    feature_tps <- (evaluate_residual_shared_information(residual_seqs[[1]], residual_seqs[[2]]))$feature_shared / nrow(residual_seqs[[1]]) * fps / log(2.0);

    outliers_and_throughputs <- cbind(outliers, feature_tps)

    return(outliers_and_throughputs)
}

plot_outliers <- function(limit = 6, fps = 120) {
    files <- length(dir("aligneddata", "^[[:digit:]]+_ali_[[:digit:]]+.txt$"))
    all_outliers <- c()

    for (i in 1:(files/2)) {
        j = 2 * i - 1
        k = j + 1

        outliers1 <- detect_outliers(j, k, limit, fps)
        outliers2 <- detect_outliers(k, j, limit, fps)

        all_outliers <- rbind(all_outliers, outliers1, outliers2)
    }

    max_outliers <- max(all_outliers[,1])
    max_tp <- floor(max(all_outliers[,2]))
    min_tp <- floor(min(all_outliers[,2]))

    par(mfrow = c(2, 2))
    outlierfit <- lm(all_outliers[,2] ~ all_outliers[,1])
    plot(all_outliers, xlim = c(0, max_outliers), ylim = c(min_tp, max_tp),
        main = "Outliers vs. feature throughput",
        xlab = "Number of outliers", ylab = "Throughput")
    abline(outlierfit, col = "red")

    print(anova(outlierfit))

    tp_means <- matrix(data = 0, nrow = max_outliers+1, ncol = 3,
        dimnames = list(c(0:max_outliers), c("mean tp", "sdev", "occurrences")))
    for (i in 0:max_outliers) {
        ivals = all_outliers[all_outliers[,1] == i,2]
        tp_means[i+1, 1] <- mean(ivals)
        tp_means[i+1, 2] <- sd(ivals)
        tp_means[i+1, 3] <- length(ivals)
    }

    barplot(tp_means[,3], names.arg = 0:max_outliers,
        main = "Number of occurrences",
        ylab = "Occurrences", xlab = "Number of outliers")

    plot(0:max_outliers, tp_means[,1], col = "red", pch = 18,
        main = "Throughput mean", xlab = "Number of outliers",
        ylab = "Throughput mean")
    plot(0:max_outliers, tp_means[,2], col = "blue", pch = 18,
        main = "Throughput standard deviation",
        xlab = "Number of outliers", ylab = "Throughput sdev")
}
