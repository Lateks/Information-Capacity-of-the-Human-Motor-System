step = 2

# Return residuals from the AR(2) model as a list where each element is the
# residual vector for the corresponding feature, given the observed sequence.
residuals <- function(sequence) {
    n <- nrow(sequence) - step
    sequence_predictors <- cbind(sequence[1:n,],sequence[1:n+step-1,])
    observed_sequence <- sequence[1:n+step,]

    frames <- nrow(observed_sequence)
    num_features <- ncol(observed_sequence)
    resid = matrix(nrow = frames, ncol = num_features)

    for (i in 1:num_features)
        resid[, i] <- get_residuals_for_feature(i, sequence_predictors,
                                                    observed_sequence)
    return(resid)
}

# Calculates the residuals for a given feature in the regression
# model.
#
# Parameters:
# - feat        the column number of the feature
# - predictors  AR model predictors
# - observed    the observed sequence
get_residuals_for_feature <- function(feat, predictors, observed) {
    frames <- nrow(observed)
    if (ncol(predictors) == ncol(observed)) # extract regressors depending on input size
        predictors <- qr(cbind(matrix(1,frames,1), predictors[, feat]))
    else
        predictors <- qr(cbind(matrix(1,frames,1), predictors[, feat],
            predictors[, ncol(observed) + feat]))
    return(qr.resid(predictors, observed[, feat]))
}

# Evaluates shared information according to residuals of residuals.
#
# Parameters:
# - residuals_a, residuals_b    residual sequences (data frames or matrices)
evaluate_residual_shared_information <- function(residuals_a, residuals_b) {
    total_shared <- 0
    total_RSS <- 0
    total_RSS_residual <- 0
    feature_shared <- array(0, ncol(residuals_a))
    n <- nrow(residuals_a)

    for (i in 1:ncol(residuals_a)) { # calculate shared information by feature
        RSS <- sum(residuals_a[, i]^2)
        total_RSS <- total_RSS + RSS

        RSS_residual <- sum_squared_residuals(residuals_a[, i], residuals_b[, i])
        total_RSS_residual <- total_RSS_residual + RSS_residual

        feature_shared[i] <- (n/2*log(RSS/RSS_residual) - log(n)/2)
        total_shared <- total_shared + feature_shared[i]
    }

    return(list(total_shared = total_shared, feature_shared = feature_shared,
        RSS = total_RSS, RSS_conditional = total_RSS_residual))
}

# Fit a linear model seq_a ~ seq_b and return the sum of squared
# residuals from the model.
sum_squared_residuals <- function(seq_a, seq_b) {
    lmfit = lm(seq_a ~ seq_b)
    residuals = resid(lmfit)
    return(sum(residuals^2))
}

# Performs PCA on data and returns a matrix of eigenvectors.
# The number of eigenvectors is chosen so that 90% of variance
# is covered.
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
