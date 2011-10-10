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

use_warnings = FALSE
twopie = 2*pi*exp(1)

# Return stochastic complexity of sequence r with side information v
evalSC <- function(v, r) {
    frames <- nrow(r)
    feature_code_length <- array(0,ncol(r))
    total_code_length <- 0
    residual_sums = list()
    for (i in 1:ncol(r)) {
        # extract regressors depending on the input size
        if (ncol(v) == ncol(r)) {
            Aplus <- qr(cbind(matrix(1,frames,1), v[, i]))
            k = 3
        }
        if (ncol(v) == 2 * ncol(r)) {
            Aplus <- qr(cbind(matrix(1,frames,1), v[, i], v[, ncol(r)+i]))
            k = 4
        }
        if (ncol(v) == 3 * ncol(r)) {
            Aplus <- qr(cbind(matrix(1,frames,1), v[, i], v[, ncol(r)+i], 
                        v[, 2*ncol(r)+i]))
            k = 5
        }
        if (ncol(v) == 4 * ncol(r))
        {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i], v[, ncol(r)+i], 
                        v[, 2*ncol(r)+i],
                        v[, 3*ncol(r)+i]))
            k = 6
        }

        res <- qr.resid(Aplus, r[, i])
        
        # warning about possibly too low variance
        if (use_warnings)
            if (sum(res^2) < .01) print(c(i,sum(res^2)))

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

# Return residuals from AR(2) as a list where each element is the
# residual vector for the corresponding feature
residuals <- function(v,r) {
    frames <- nrow(r)
    residuals = list()
    # extract regressors depending on the input size
    for (i in 1:ncol(r)) {
        if (ncol(v) == ncol(r))
            Aplus <- qr(cbind(matrix(1,frames,1), v[, i]))
        if (ncol(v) == 2 * ncol(r))
            Aplus <- qr(cbind(matrix(1,frames,1), v[, i], v[, ncol(r)+i]))

        residuals[[i]] <- qr.resid(Aplus, r[, i])  
    }
    return(residuals)
}

shared_information_by_residuals <- function(a,b,index) {
    n <- min(nrow(a), nrow(b))-index
    total_shared <- 0
    feature_shared <- array(0,ncol(a))
    
    resa <- residuals(cbind(a[1:n,],a[1:n+index-1,]), a[1:n+index,])
    resb <- residuals(cbind(b[1:n,],b[1:n+index-1,]), b[1:n+index,])

    for (i in 1:length(resa)) {
        lm.ra = lm(resa[[i]] ~ resb[[i]])
        resaf = resid(lm.ra)

        rssa <- sum(resa[[i]]^2)
        rssaf <- sum(resaf^2)
        
        shared_information <- (n/2*log(rssa/rssaf) - log(n)/2)
        
        feature_shared[i] <- shared_information
        total_shared <- total_shared + shared_information
    }

    return(list(n = n, total_shared = total_shared, feature_shared = feature_shared))
}

# Complexity of a single (multivariate) sequence
SC <- function(a) {
    n <- nrow(a)-2
    if (n < 1)
        return(list(n = 0,k = 1,cl = 0,rate = 0))
    
    SC <- evalSC(cbind(a[1:n,],a[1:n+1,]), a[1:n+2,])
    complexity <- SC$total
    feature_complexity <- SC$feature
    residual_sums <- SC$residual_sums
    return(list(n = n, complexity = complexity, rate = complexity/n,
        feature_complexity = feature_complexity, residual_sums = residual_sums))
}

# Complexity of a (multivariate) sequence given another
SCcond <- function(a,b) {
    n <- min(nrow(a), nrow(b))-2
    SC <- evalSC(cbind(a[1:n,],a[1:n+1,],b[1:n+2,]), a[1:n+2,])
    complexity <- SC$total
    feature_complexity <- SC$feature
    residual_sums <- SC$residual_sums
    return(list(n = n, complexity = complexity, rate = complexity/n,
        feature_complexity = feature_complexity, residual_sums = residual_sums))
}

# Removes duplicated rows from both sequences.
remove_duplicate_frames <- function(a, b, method) {
    skipa <- matrix(FALSE, nrow(a))
    skipa <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) < 0.001
    if (method == 2) {
        skipb <- matrix(FALSE, nrow(b))
        skipb <- rowSums((b[2:nrow(b),]-b[1:(nrow(b)-1),])^2) < 0.001
        skip <- skipa | skipb
    }
    if (method == 1) {
        skip <- skipa
    }
    a <- a[!skip,]
    b <- b[!skip,]

    return(list(a, b))
}

normalize_features <- function(a) {
    a <- t(apply(a,1,'-',apply(a,2,mean)))
    a <- t(apply(a,1,'/',sqrt(apply(a,2,var))))
    a[is.nan(a)]<-0

    return(a)
}

# Add noise to the features of a matrix.    
add_noise_to_features <- function(a, noise_coeff) {
    features <- ncol(a)
    for (k in 1:features) {
        feature_sdev <- sqrt(var(a[,k]))
        a[,k] <- a[,k] + noise_coeff * feature_sdev * rnorm(nrow(a))
    }
    return(a)
}

# Returns a matrix containing the chosen most significant eigenvectors
pca <- function(data) {
    principal_components <- prcomp(data, retx = TRUE, center = TRUE, scale. = TRUE)

    sum <- 0
    eigenvectors <- 0
    threshold <- 0.9*sum(principal_components$sdev**2)
    for (k in 1:length(principal_components$sdev)) {
        sum <- sum + principal_components$sdev[k]**2
        eigenvectors <- k
        if (sum >= threshold)
            break
    }

    eigenvectors <- principal_components$rotation[,1:eigenvectors]

    return(eigenvectors)
}

# Returns throughput determined by residuals of residuals.
residual_throughput <- function(a, b, fps = 120, features = FALSE, index = 2) {
    lena <- nrow(a)

    shared_information <- shared_information_by_residuals(a, b, index)
    throughput <- shared_information$total_shared/lena * fps/log(2.0)

    if (features) {
        feature_throughputs <- shared_information$feature_shared/lena * fps/log(2.0)
        return(list(throughput = throughput,
            feature_throughputs = feature_throughputs))
    }

    return(list(throughput = throughput))
}

# Returns throughput calculated in the manner described in the original
# arxiv paper draft.
original_throughput <- function(a, b, fps = 120, features = FALSE) {
    lena <- nrow(a)
    SC_a <- SC(a)
    SC_a_b <- SCcond(a, b)

    complexity_of_a <- SC_a$complexity
    feature_complexity_of_a <- SC_a$feature_complexity

    complexity_of_a_cond_b <- SC_a_b$complexity
    feature_complexity_of_a_cond_b <- SC_a_b$feature_complexity

    throughput <- (complexity_of_a - complexity_of_a_cond_b)/lena*fps/log(2.0)
    
    if (features) {
        feature_throughputs <- array(0,ncol(a))
        for (k in 1:ncol(a))
            feature_throughputs[k] <- (feature_complexity_of_a[k] -
                feature_complexity_of_a_cond_b[k])/lena*fps/log(2.0)

        return(list(throughput = throughput,
            feature_throughputs = feature_throughputs))
    }

    return(list(throughput = throughput))
}

# Calculate throughput for a given pair of matrices.
throughput <- function(a, b, fps = 120, pca = FALSE, residuals = FALSE,
    features = FALSE, index = 2) {

    lena <- nrow(a)

    if (pca) {
        eigenvectors <- pca(a)
        a <- a %*% eigenvectors
        b <- b %*% eigenvectors
    }
    
    if (residuals)
        return(residual_throughput(a, b, fps = fps, features = features, index = index))

    return(original_throughput(a, b, fps = fps , features = features))
}

# Calculate throughput for a given pair of files.
#
# See the dir_throughput function for descriptions of optional parameters.
pair_throughput <- function(filename1, filename2, fps = 120, pca = FALSE,
    amc = FALSE, residuals = FALSE, features = FALSE, warnings = FALSE,
    noise = 0, index = 2, method = 1) {
    
    use_warnings <<- warnings

    a <- read.table(filename1)
    b <- read.table(filename2)

    # Eliminate features with zero variance in AMC format data, usually the
    # fingers of both hands.
    if (amc) {
        a$V34=NULL; a$V46=NULL;
        b$V34=NULL; b$V46=NULL;
    }

    data <- remove_duplicate_frames(a, b, method)
    a <- data[[1]]
    b <- data[[2]]

    if (noise > 0) {
        a <- add_noise_to_features(a, noise)
        b <- add_noise_to_features(b, noise)
    }

    a <- normalize_features(a)
    b <- normalize_features(b)

    use_warnings <<- FALSE
    return(throughput(a, b, fps = fps, pca = pca, residuals = residuals,
        features = features, index = index))
}

# Calculate throughput for all aligned sequence pairs in the directory (by
# default the current working directory).
#
# The function assumes that the original amc or coordinate files are in the
# working directory. (Note: all files should be for the same test subject!)
# Files should be named <sequence number>.txt, e.g. "15.txt".
#
# A subdirectory named 'aligneddata' should contain the aligned data
# in files that are named in the following way:
#         <seq1 #>_ali_<seq2 #>.txt
# Example: "14_ali_15.txt".
#
# Parameters:
# - fps        frames per second
# - pca        use principal components analysis
# - amc        input data is in AMC format
# - residuals  use residuals to determine complexities
# - features   also calculate throughput for individual features
# - warnings   print warnings of low residual variance (to debug NaN results)
# - noise      noise coefficient (default 0)
# - index      the index to customize AR(2) of "residual of residuals" method
# - method     the method to remove duplicates (default 1)
dir_throughput <- function(fps = 120, pca = FALSE, amc = FALSE, residuals = FALSE,
    features = FALSE, warnings = FALSE, noise = 0, index = 2, method = 1) {

    if (amc) {
        files <- dir(".", "^[[:digit:]]+.amc$")
    } else
        files <- dir(".", "^[[:digit:]]+.txt$")
    sequences <- length(files)
    
    total_throughputs <- array(0, c(sequences/2, 2))
    feature_throughputs <- array(list(NULL), c(sequences/2, 2))
	
    for (i in 1:(sequences/2)) {
        j = 2 * i - 1
        k = j + 1
        file1 <- sprintf("aligneddata/%d_ali_%d.txt", j, k)
        file2 <- sprintf("aligneddata/%d_ali_%d.txt", k, j)

        results <- pair_throughput(file1, file2, fps = fps, pca = pca,
            amc = amc, res = residuals, features = features, warnings = warnings,
            noise = noise, index = index, method = method)
        inverse_results <- pair_throughput(file2, file1, fps = fps, pca = pca,
            amc = amc, res = residuals, features = features, warnings = warnings,
            noise = noise, index = index, method = method)

        total_throughputs[i,1] <- results$throughput
        total_throughputs[i,2] <- inverse_results$throughput

        if (features) {
            feature_throughputs[i,1] <- list(results$feature_throughputs)
            feature_throughputs[i,2] <- list(inverse_results$feature_throughputs)
        }
    }
    
    if (features) {
        return(list(total = total_throughputs, feature = feature_throughputs))
    } else
        return(list(total = total_throughputs))
}