step = 2

# Given a result list, FPS count and sequence length, generates a
# vector with throughput, shared information in bits, total sum of
# residuals, total conditional sum of residuals and the residual
# quotient. (Vector length is 5.)
#
# Parameters:
# - results             a result list given by evaluate_residual_shared_information
# - fps                 frames per second (to calculate throughut)
# - seq_length          sequence length (after alignment and duplicate removal)
construct_result_vector <- function(results, fps, seq_length) {
    quotient <- results$RSS / results$RSS_conditional
    throughput <- results$total_shared / seq_length * fps / log(2.0)
    shared_information_bits <- results$total_shared / log(2.0)

    return(c(throughput, shared_information_bits, results$RSS,
        results$RSS_conditional, quotient))
}

# Return residuals from the AR(2) model as a list where each element is the
# residual vector for the corresponding feature, given the observed sequence
# and some predictor features (shifted data from the observed sequence).
evaluate_residuals <- function(sequence_predictors, observed_sequence) {
    frames <- nrow(observed_sequence)
    num_features <- ncol(observed_sequence)
    residuals = matrix(nrow = frames, ncol = num_features)
    # extract regressors depending on the input size
    for (i in 1:num_features) {
        if (ncol(sequence_predictors) == num_features)
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i]))
        if (ncol(sequence_predictors) == 2 * num_features)
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i],
                sequence_predictors[, num_features+i]))

        residuals[, i] <- qr.resid(predictors, observed_sequence[, i])
    }
    return(residuals)
}

# Evaluates shared information using the "residuals of residuals" method.
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

# Returns shared information calculated with the "residuals of residuals" method.
# Parameters:
# - a, b    sequences
shared_information_by_residuals <- function(a,b) {
    n <- nrow(a) - step

    residuals_a <- evaluate_residuals(cbind(a[1:n,],a[1:n+step-1,]), a[1:n+step,])
    residuals_b <- evaluate_residuals(cbind(b[1:n,],b[1:n+step-1,]), b[1:n+step,])

    results <- evaluate_residual_shared_information(residuals_a, residuals_b)

    return(list(n = n, total_shared = results$total_shared,
        feature_shared = results$feature_shared, RSS = results$RSS,
        RSS_conditional = results$RSS_conditional))
}

# Removes duplicated rows from both sequences a and b
# either symmetrically or asymmetrically.
#
# Returns the altered sequences as a list with two elements.
#
# Parameters:
# - a, b        sequences (order matters if asymmetric removal
#               is used)
remove_duplicate_frames <- function(a, b) {
    skipa <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) == 0
    skipa <- c(FALSE, skipa)
    skip <- skipa

    a <- a[!skip,]
    b <- b[!skip,]

    return(list(a, b))
}

# Performs normalization on the features of sequence a.
# Returns the altered sequence.
normalize_features <- function(a) {
    a <- t(apply(a,1,'-',apply(a,2,mean)))
    a <- t(apply(a,1,'/',sqrt(apply(a,2,var))))
    a[is.nan(a)]<-0

    return(a)
}

choose_number_of_eigenvectors <- function(sdevs) {
    sum <- 0
    eigenvectors <- 0
    threshold <- 0.9*sum(sdevs**2)

    for (k in 1:length(sdevs)) {
        sum <- sum + sdevs[k]**2
        eigenvectors <- k
        if (sum >= threshold)
            break
    }

    return(eigenvectors)
}

# Returns a list containing both of the given matrices after
# being multiplied with the chosen eigenvectors (calculated
# from the matrix data).
#
# Parameters:
# - data        data matrix to perform PCA to
# - side_data   another matrix to multiply with the same eigenvectors
pca <- function(data, side_data) {
    # Normalization is done here instead of passing center = TRUE and
    # scale. = TRUE to prcomp because multiplication with eigenvectors
    # is done "manually" (outside of prcomp) and therefore at least mean
    # normalization has to be done on the original data anyway.
    data <- normalize_features(data)
    side_data <- normalize_features(side_data)

    principal_components <- prcomp(data)
    num_eigenvecs <- choose_number_of_eigenvectors(principal_components$sdev)

    eigenvectors <- principal_components$rotation[,1:num_eigenvecs]
    reduced_data <- data %*% eigenvectors
    reduced_side_data <- side_data %*% eigenvectors

    return(list(reduced_data, reduced_side_data))
}
