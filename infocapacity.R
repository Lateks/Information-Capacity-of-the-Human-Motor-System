#####################################################################
#
# Original:
# R code for reproducing the experiments in 
# (Roos, Oulasvirta, "An Extended Framework for Measuring the
# Information Capacity of the Human Motor System")
#
# code by Teemu Roos, Feb 15, 2011
# teemu.roos (at) cs.helsinki.fi
#
# Current version:
# Modified for more general use and added PCA functionality
# by Laura Leppänen (research assistant), September 2011
# laura.leppanen (at) cs.helsinki.fi
#
# New method of calculating throughput using residuals of residuals
# added by Arttu Modig (research assistant), September 2011
# arttu.modig (at) aalto.fi
#
# Data needs to be aligned temporally before being evaluated
# for information capacity. To align, we recommend Canonical
# Time Warping, from www.humansensing.cs.cmu.edu/projects/ctwCode.html
#
#####################################################################

options(width=Sys.getenv("COLUMNS"))

twopie = 2*pi*exp(1)

# Given a result list, FPS count and sequence length, generates a
# vector with throughput, shared information in bits, total sum of
# residuals, total conditional sum of residuals and the residual
# quotient. (Vector length is 5.)
#
# Parameters:
# - results             a result list given by evaluate_residual_shared_information
# - fps                 frames per second (to calculate throughut)
# - seq_length          sequence length (after alignment and duplicate removal)
construct_result_vector <- function(results, fps, seq_length)
{
    quotient <- results$RSS / results$RSS_conditional
    throughput <- results$total_shared / seq_length * fps / log(2.0)
    shared_information_bits <- results$total_shared / log(2.0)

    return(c(throughput, shared_information_bits, results$RSS, results$RSS_conditional, quotient))
}

# Return stochastic complexity of the observed sequence given predictor features
# (sometimes containing side information from another sequence, but usually only
# shifted data from the same sequence).
evaluate_SC <- function(sequence_predictors, observed_sequence) {
    frames <- nrow(observed_sequence)
    num_features <- ncol(observed_sequence)
    feature_code_length <- array(0,ncol(observed_sequence))
    total_code_length <- 0
    residual_sums = list()
    for (i in 1:num_features) {
        # extract regressors depending on the input size
        if (ncol(sequence_predictors) == num_features) {
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i]))
            k = 3
        }
        if (ncol(sequence_predictors) == 2 * num_features) {
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i], sequence_predictors[, num_features+i]))
            k = 4
        }
        if (ncol(sequence_predictors) == 3 * num_features) {
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i], sequence_predictors[, num_features+i],
                        sequence_predictors[, 2*num_features+i]))
            k = 5
        }
        if (ncol(sequence_predictors) == 4 * num_features) {
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i], sequence_predictors[, num_features+i],
                        sequence_predictors[, 2*num_features+i],
                        sequence_predictors[, 3*num_features+i]))
            k = 6
        }

        res <- qr.resid(predictors, observed_sequence[, i])

        residual_sum <- sum(res^2)
        residual_sums[[i]] <- residual_sum

        # Rissanen's classic two-part MDL code-length (feel free
        # to replace by more advanced universal code).
        MDL <-  frames/2*log(twopie*residual_sum/frames) + k/2*log(frames)

        feature_code_length[i] <- MDL
        total_code_length <- total_code_length + MDL
    }
    return(list(total = total_code_length, feature = feature_code_length,
        residual_sums = residual_sums))
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
            predictors <- qr(cbind(matrix(1,frames,1), sequence_predictors[, i], sequence_predictors[, num_features+i]))

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
# - index   desired size of the second step in the AR model
shared_information_by_residuals <- function(a,b,index) {
    n <- nrow(a) - index

    residuals_a <- evaluate_residuals(cbind(a[1:n,],a[1:n+index-1,]), a[1:n+index,])
    residuals_b <- evaluate_residuals(cbind(b[1:n,],b[1:n+index-1,]), b[1:n+index,])

    results <- evaluate_residual_shared_information(residuals_a, residuals_b)

    return(list(n = n, total_shared = results$total_shared,
        feature_shared = results$feature_shared, RSS = results$RSS,
        RSS_conditional = results$RSS_conditional))
}

# Complexity of a single (multivariate) sequence
# Parameters:
# - a       sequence
# - index   desired size of the second step in the AR model
SC <- function(a, index) {
    n <- nrow(a)-index
    if (n < 1)
        return(list(n = 0,k = 1,cl = 0,rate = 0))

    SC <- evaluate_SC(cbind(a[1:n,],a[1:n+index-1,]), a[1:n+index,])
    complexity <- SC$total
    feature_complexity <- SC$feature
    residual_sums <- SC$residual_sums
    return(list(n = n, complexity = complexity, rate = complexity/n,
        feature_complexity = feature_complexity, residual_sums = residual_sums))
}

# Complexity of a (multivariate) sequence given another
# Parameters:
# - a       primary sequence
# - b       side information equence
# - index   desired size of the second step in the AR model
SCcond <- function(a,b, index) {
    n <- min(nrow(a), nrow(b))-index
    SC <- evaluate_SC(cbind(a[1:n,],a[1:n+index-1,],b[1:n+index,]), a[1:n+index,])
    complexity <- SC$total
    feature_complexity <- SC$feature
    residual_sums <- SC$residual_sums
    return(list(n = n, complexity = complexity, rate = complexity/n,
        feature_complexity = feature_complexity, residual_sums = residual_sums))
}

# Removes duplicated rows from both sequences a and b
# either symmetrically or asymmetrically.
#
# Returns the altered sequences as a list with two elements.
#
# Parameters:
# - a, b        sequences (order matters if asymmetric removal
#               is used)
# - symmetric   boolean value indicating removal type
#               (TRUE = symmetric, FALSE = asymmetric)
remove_duplicate_frames <- function(a, b, symmetric = FALSE) {
    skipa <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) == 0
    skipa <- c(FALSE, skipa)
    if (symmetric) {
        skipb <- matrix(FALSE, nrow(b))
        skipb <- rowSums((b[2:nrow(b),]-b[1:(nrow(b)-1),])^2) == 0
        skipb <- c(FALSE, skipb)
        skip <- skipa | skipb
    }
    else
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

# Add noise to the features of a matrix.
# Parameters:
# - a               seuqence matrix
# - noise_coeff     coefficient for standard deviation
add_noise_to_features <- function(a, noise_coeff) {
    features <- ncol(a)
    for (k in 1:features) {
        feature_sdev <- sqrt(var(a[,k]))
        a[,k] <- a[,k] + noise_coeff * feature_sdev * rnorm(nrow(a))
    }
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

# Returns throughput determined by residuals of residuals.
#
# Parameters:
# - a, b        sequences
# - fps         frames per second
# - features    return throughput and shared information for individual
#               features
# - index       desired size of the second step in the AR model
residual_throughput <- function(a, b, fps = 120, features = FALSE, index = 2) {
    lena <- nrow(a)

    shared <- shared_information_by_residuals(a, b, index)
    total_shared <- shared$total_shared

    if (features) {
        feature_shared <- shared$feature_shared
        feature_throughputs <- feature_shared / lena * fps / log(2.0)
        return(list(total_shared = total_shared, feature_shared = feature_shared,
            feature_throughputs = feature_throughputs,
            RSS = shared$RSS, RSS_conditional = shared$RSS_conditional))
    }

    return(list(total_shared = total_shared, RSS = shared$RSS,
        RSS_conditional = shared$RSS_conditional))
}

# Returns throughput calculated in the manner described in the original
# arxiv paper draft.
#
# Parameters:
# - a, b        sequences
# - fps         frames per second
# - features    return throughput and shared information for individual
#               features
# - index       desired size of the second step in the AR model
original_throughput <- function(a, b, fps = 120, features = FALSE, index = 2) {
    SC_a <- SC(a, index)
    SC_a_b <- SCcond(a, b, index)

    complexity_of_a <- SC_a$complexity
    feature_complexity_of_a <- SC_a$feature_complexity

    complexity_of_a_cond_b <- SC_a_b$complexity
    feature_complexity_of_a_cond_b <- SC_a_b$feature_complexity

    shared_information <- complexity_of_a - complexity_of_a_cond_b

    RSS_a <- sum(unlist(SC_a$residual_sums))
    RSS_a_cond_b <- sum(unlist(SC_a_b$residual_sums))

    if (features) {
        length_a <- nrow(a)
        feature_shared_informations <- array(0,ncol(a))
        feature_throughputs <- array(0,ncol(a))
        for (k in 1:ncol(a)) {
            feature_shared_informations[k] <- (feature_complexity_of_a[k] -
                feature_complexity_of_a_cond_b[k])
            feature_throughputs[k] <- (feature_complexity_of_a[k] -
                feature_complexity_of_a_cond_b[k]) / length_a * fps / log(2.0)
        }

        return(list(total_shared = shared_information,
            feature_shared = feature_shared_informations,
            feature_throughputs = feature_throughputs, RSS = RSS_a,
            RSS_conditional = RSS_a_cond_b))
    }

    return(list(total_shared = shared_information, RSS = RSS_a,
        RSS_conditional = RSS_a_cond_b))
}

# Calculate throughput for a given pair of matrices.
#
# Parameters:
# - a, b        sequences
# - fps         frames per second
# - residuals   use the "residuals of residuals" method
# - features    return throughput and shared information for individual
#               features
# - index       desired size of the second step in the AR model
throughput <- function(a, b, fps = 120, pca = FALSE, residuals = TRUE,
    features = FALSE, index = 2) {

    if (pca) {
        reduced <- pca(a, b)
        a <- reduced[[1]]
        b <- reduced[[2]]
    }

    if (residuals)
        return(residual_throughput(a, b, fps = fps,
            features = features, index = index))

    return(original_throughput(a, b, fps = fps , features = features, index = index))
}

# Loads the aligned data file "seqnum1_ali_seqnum2.txt" as a data frame
# and returns it.
load_aligned <- function(seqnum1, seqnum2) {
    return(read.table(sprintf("aligneddata/%d_ali_%d.txt", seqnum1, seqnum2)))
}

# Calculate throughput for a given pair of files.
#
# See the dir_throughput function for descriptions of optional parameters.
pair_throughput <- function(seqnum1, seqnum2, fps = 120, pca = FALSE,
    residuals = TRUE, features = FALSE, noise = 0, index = 2,
    symmetric = FALSE) {

    a <- load_aligned(seqnum1, seqnum2)
    b <- load_aligned(seqnum2, seqnum1)

    data <- remove_duplicate_frames(a, b, symmetric)
    a <- data[[1]]
    b <- data[[2]]

    if (noise > 0)
        a <- add_noise_to_features(a, noise)

    results <- throughput(a, b, fps = fps, pca = pca, residuals = residuals,
        features = features, index = index)
    if (features)
        return(list(results = construct_result_vector(results, fps, nrow(a)),
            feature_throughputs = results$feature_throughputs))
    return(list(results = construct_result_vector(results, fps, nrow(a))))
}

# For a given pair of files, calculate residual sums and throughputs with
# different indices for the second AR time step.
#
# See dir_throughput for optional parameters. The index parameter that
# supplies the step size for the AR model is replaced with three new
# parameters:
#
# start_step    the starting value for the AR model step size (default 2)
# end_step      the end value for the AR model step size (default 120)
# interval      the interval between step sizes in the trials (default 1)
evaluate_step_series <- function(filename1, filename2, fps = 120, pca = FALSE,
    residuals = TRUE, features = FALSE, noise = 0,
    symmetric = FALSE, start_step = 2, end_step = 120, interval = 1) {

    size <- (end_step-start_step)/interval+1
    results <- array(0, c(size, 5));
    colnames(results) <- c("index","RSS","RSS_conditional","RSS / RSS_cond","TP")
    i <- 1

    for (index in seq(start_index, end_index, interval)) {
        print(index)

        TP <- pair_throughput(filename1, filename2, fps = fps, pca = pca,
            residuals = residuals, features = features,
            noise = noise, index = index, symmetric = symmetric)

        results[i,] <- TP$results
        i <- i+1
    }

    return(results)
}

# Calculate throughput for all aligned sequence pairs in the directory (by
# default the current working directory).
#
# The function assumes that the original coordinate files are in the
# working directory. (Note: all files should be for the same test subject!)
# Files should be named <sequence number>.txt, e.g. "15.txt".
#
# The aligned data should reside in the subdirectory "aligneddata"
# in files that are named in the following way:
#         <seq1 #>_ali_<seq2 #>.txt
# Example: "14_ali_15.txt".
#
# Parameters:
# - fps        frames per second
# - pca        use principal components analysis
# - residuals  use residuals of residuals to determine complexities
# - features   also calculate throughput for individual features
# - noise      noise coefficient (default 0)
# - index      the step variable for the AR(2) model in the "residuals of residuals"
#              method
# - symmetric  remove duplicates symmetrically (default is asymmetric)
dir_throughput <- function(fps = 120, pca = FALSE, residuals = TRUE,
    features = FALSE, noise = 0, index = 2,
    symmetric = FALSE) {

    files <- dir(".", "^[[:digit:]]+.txt$")
    sequences <- length(files)

    col_names <- c("TP", "shared", "RSS", "RSS_cond", "quotient")
    total_throughputs <- matrix(nrow = sequences, ncol = 5,
        dimnames = list(1:sequences, col_names))
    rownames = c()

    feature_throughputs <- array(list(NULL), c(sequences/2, 2))

    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        rownames <- c(rownames, sprintf("(%d, %d)", j, k))
        rownames <- c(rownames, sprintf("(%d, %d)", k, j))

        results <- pair_throughput(j, k, fps = fps, pca = pca,
            residuals = residuals, features = features,
            noise = noise, index = index, symmetric = symmetric)
        inverse_results <- pair_throughput(k, j, fps = fps, pca = pca,
            residuals = residuals, features = features,
            noise = noise, index = index, symmetric = symmetric)

        total_throughputs[j,] <- results$results
        total_throughputs[k,] <- inverse_results$results

        if (features) {
            feature_throughputs[i,1] <- list(results$feature_throughputs)
            feature_throughputs[i,2] <- list(inverse_results$feature_throughputs)
        }
    }
    rownames(total_throughputs) <- rownames

    if (features) {
        return(list(total = total_throughputs, feature = feature_throughputs))
    } else
        return(list(total = total_throughputs))
}
