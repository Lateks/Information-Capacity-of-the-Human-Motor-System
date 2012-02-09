step = 2

# Returns shared information calculated with the "residuals of residuals" method.
# Parameters:
# - a, b    sequences
shared_information_by_residuals <- function(a,b) {
    residuals_a <- residuals(a)
    residuals_b <- residuals(b)

    results <- evaluate_residual_shared_information(residuals_a, residuals_b)

    return(append(results, n = nrow(a) - step))
}

# Return residuals from the AR(2) model as a list where each element is the
# residual vector for the corresponding feature, given the observed sequence
# and some predictor features (shifted data from the observed sequence).
residuals <- function(sequence) {
    n <- nrow(sequence) - step
    sequence_predictors <- cbind(sequence[1:n,],sequence[1:n+step-1,])
    observed_sequence <- sequence[1:n+step,]

    frames <- nrow(observed_sequence)
    num_features <- ncol(observed_sequence)
    residuals = matrix(nrow = frames, ncol = num_features)

    for (i in 1:num_features) { # extract regressors depending on the input size
        if (ncol(sequence_predictors) == num_features)
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i]))
        if (ncol(sequence_predictors) == 2 * num_features)
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i],
                sequence_predictors[, num_features + i]))
        residuals[, i] <- qr.resid(predictors, observed_sequence[, i])
    }
    return(residuals)
}

# Evaluates shared information using the "residuals of residuals" method.
#
# Parameters:
# - residuals_a, residuals_b    residual sequences (data frames or matrices)
evaluate_residual_shared_information <- function(residuals_a, residuals_b) {
    total_shared <- 0
    total_RSS <- 0
    total_RSS_residual <- 0
    feature_shared <- array(0, ncol(residuals_a))
    n <- nrow(residuals_a)

    for (i in 1:ncol(residuals_a)) {
        lm.ra = lm(residuals_a[, i] ~ residuals_b[, i])
        residuals = resid(lm.ra)

        RSS <- sum(residuals_a[, i]^2)
        RSS_residual <- sum(residuals^2)

        total_RSS <- total_RSS + RSS
        total_RSS_residual <- total_RSS_residual + RSS_residual

        shared_information <- (n/2*log(RSS/RSS_residual) - log(n)/2)

        feature_shared[i] <- shared_information
        total_shared <- total_shared + shared_information
    }

    return(list(total_shared = total_shared, feature_shared = feature_shared,
        RSS = total_RSS, RSS_conditional = total_RSS_residual))
}

# Removes duplicated rows from both sequences a and b.
#
# Returns the altered sequences as a list with two elements.
remove_duplicate_frames <- function(a, b) {
    skip <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) == 0
    skip <- c(FALSE, skip)

    a <- a[!skip,]
    b <- b[!skip,]

    return(list(a, b))
}

# Performs PCA on data and returns a matrix of eigenvectors.
pca <- function(data) {
    data <- normalize_features(data)

    principal_components <- prcomp(data)
    num_eigenvecs <- choose_number_of_eigenvectors(principal_components$sdev)

    eigenvectors <- principal_components$rotation[,1:num_eigenvecs]
    return(eigenvectors)
}

# Performs normalization on the features of sequence a.
# Returns the altered sequence.
normalize_features <- function(a) {
    a <- t(apply(a,1,'-',apply(a,2,mean)))
    a <- t(apply(a,1,'/',sqrt(apply(a,2,var))))
    a[is.nan(a)]<-0

    return(a)
}

# Given a list of standard deviations (from the prcomp function
# output), choose a number of eigenvectors in such a way that
# 90% of variance is covered.
choose_number_of_eigenvectors <- function(sdevs) {
    threshold <- 0.9*sum(sdevs**2)
    sum <- 0
    for (k in 1:length(sdevs)) {
        sum <- sum + sdevs[k]**2
        if (sum >= threshold)
            return(k)
    }
}
