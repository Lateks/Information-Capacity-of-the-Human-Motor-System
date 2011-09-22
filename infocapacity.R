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
# by Laura Leppänen (research assistant), Sep 9, 2011
# laura.leppanen (at) cs.helsinki.fi
#
# Data needs to be aligned temporally before being evaluated
# for information capacity. To align, we recommend Canonical
# Time Warping, from www.humansensing.cs.cmu.edu/projects/ctwCode.html
#
#####################################################################

options(width=Sys.getenv("COLUMNS"));
twopie = 2*pi*exp(1);

# Return stochastic complexity of sequence r with side information v
evalSC <- function(v, r, war = FALSE) {
    n <- nrow(r); # number of frames
    feature_code_length <- array(0,ncol(r));
    total_code_length <- 0;
    for (i in 1:ncol(r)) {
        # extract regressors depending on the input size
        if (ncol(v) == ncol(r)) {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i]));
            k = 3;
        }
        if (ncol(v) == 2 * ncol(r)) {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i], v[, ncol(r)+i]));
            k = 4;
        }
        if (ncol(v) == 3 * ncol(r)) {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i], v[, ncol(r)+i], 
                        v[, 2*ncol(r)+i]));
            k = 5;
        }
        if (ncol(v) == 4 * ncol(r))
        {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i], v[, ncol(r)+i], 
                        v[, 2*ncol(r)+i],
                        v[, 3*ncol(r)+i]));
            k = 6;
        }

        res <- qr.resid(Aplus, r[, i]);
        
        # warning about possibly too low variance
        if (war)
            if (sum(res^2) < .01) print(c(i,sum(res^2)))

        rss <- sum(res^2)

        # Rissanen's classic two-part MDL code-length (feel free
        # to replace by more advanced universal code).
        MDL <-  n/2*log(twopie*rss/n) + k/2*log(n);
        
        feature_code_length[i] <- MDL;
        total_code_length <- total_code_length + MDL;
    }
    return(list(totCL=total_code_length, featCL=feature_code_length, rss=rss));
}

# Return residuals from AR(2) as a list where each element is the
# residual vector for the corresponding feature
evalR <- function(v,r) {
    n <- nrow(r); # number of frames
    res = list();
    # extract regressors depending on the input size
    for (i in 1:ncol(r)) {
        if (ncol(v) == ncol(r)) {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i]));
            k <- 3;
        }
        if (ncol(v) == 2 * ncol(r)) {
            Aplus <- qr(cbind(matrix(1,n,1), v[, i], v[, ncol(r)+i]));
            k <- 4;
        }
        res[[i]] <- qr.resid(Aplus, r[, i]);  
    }
    return(res);
}

# Shared information of sequences a and b using complexity
# defined by residuals
shared_information <- function(a,b) {
    n <- min(nrow(a), nrow(b))-2;
    total_shared <- 0;
    feature_shared <- array(0,ncol(a));
    
    resa <- evalR(cbind(a[1:n,],a[1:n+1,]), a[1:n+2,]);
    resb <- evalR(cbind(b[1:n,],b[1:n+1,]), b[1:n+2,]);

    for (i in 1:length(resa)) {
        lm.ra = lm(resa[[i]] ~ resb[[i]]);
        resaf = resid(lm.ra);

        rssa <- sum(resa[[i]]^2);
        rssaf <- sum(resaf^2);
        
        shared_information <- (n/2*log(rssa/rssaf) - log(n)/2);
        
        feature_shared[i] <- shared_information;
        total_shared <- total_shared + shared_information;
    }

    return(list(n=n, total_shared=total_shared, feature_shared=feature_shared));
}

# Complexity of a single (multivariate) sequence
SC <- function(a, war=FALSE) {
    n <- nrow(a)-2;
    if (n < 1)
        return(list(n = 0,k = 1,cl = 0,rate = 0));
    
    SC <- evalSC(cbind(a[1:n,],a[1:n+1,]), a[1:n+2,], war=war);
    cl <- SC$totCL;
    fcl <- SC$featCL;
    rss <- SC$rss;
    return(list(n=n, cl=cl, rate=cl/n, fcl=fcl, rss=rss));
}

# Complexity of a (multivariate) sequence given another
SCcond <- function(a,b, war=FALSE) {
    n <- min(nrow(a), nrow(b))-2;
    SC <- evalSC(cbind(a[1:n,],a[1:n+1,],b[1:n+2,]), a[1:n+2,], war=war);
    cl <- SC$totCL;
    fcl <- SC$featCL;
    rss <- SC$rss;
    return(list(n=n, cl=cl, rate=cl/n, fcl=fcl, rss=rss));
}

# Removes duplicated rows from both sequences.
remove_duplicate_frames <- function(a, b) {
    skipa <- matrix(FALSE, nrow(a));
    skipa <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) < 0.001;
    skipb <- matrix(FALSE, nrow(b));
    skipb <- rowSums((b[2:nrow(b),]-b[1:(nrow(b)-1),])^2) < 0.001;
    skip <- skipa | skipb;
    a <- a[!skip,];
    b <- b[!skip,];

    return(list(a, b));
}

# Perform normalization on the features of a matrix.
normalize_features <- function(a) {
    a <- t(apply(a,1,'-',apply(a,2,mean)));
    a <- t(apply(a,1,'/',sqrt(apply(a,2,var))));
    a[is.nan(a)]<-0;

    return(a);
}

# Add noise on the features of a matrix.           
add_noise <- function(a, x) {
    nfeat <- ncol(a);
    for (k in 1:nfeat) {
        feat_var <- var(a[,k]);
        a[,k] <- a[,k] + x * feat_var * rnorm(nrow(a));
    }
    return(a);
}

# Calculate throughput for a given pair of matrices.
throughput <- function(a, b, fps = 120, pca = FALSE, res = FALSE, ftp = FALSE, war = FALSE) {
    lena <- nrow(a);

    if (pca) {
        principal_components <- prcomp(a, retx = TRUE, center = TRUE, scale. = TRUE);

        sum <- 0;
        number_of_eigenvectors <- 0;
        threshold <- 0.9*sum(principal_components$sdev**2);
        for (k in 1:length(principal_components$sdev)) {
            sum <- sum + principal_components$sdev[k]**2;
            number_of_eigenvectors <- k;
            if (sum >= threshold)
                break;
        }
        eigenvectors <- pcs$rotation[,1:number_of_eigenvectors];
        a <- a %*% eigenvectors;
        b <- b %*% eigenvectors;
    }
    
    if (res) { # Return shared information determined by residuals
        shared_information <- shared_information(a, b);
        throughput <- shared_information$total_shared/lena*120/log(2.0);
        feature_throughput <- shared_information$feature_shared/lena*120/log(2.0);
        return(list(tp=throughput, featTP=feature_throughput));
    }

    SC_a <- SC(a, war);
    SC_a_b <- SCcond(a, b, war);
    complexity_of_a <- SC_a$cl;
    feature_complexity_of_a <- SC_a$fcl;
    complexity_of_a_cond_b <- SC_a_b$cl;
    feature_complexity_of_a_cond_b <- SC_a_b$cl;

    throughput <- (complexity_of_a - complexity_of_a_cond_b)/lena*fps/log(2.0);
    
    if (ftp) {
        feature_throughputs <- array(0,ncol(a));
        for (k in 1:ncol(a)) {
            feature_throughputs[k] <- (feature_complexity_of_a[k] -
                feature_complexity_of_a_cond_b[k])/lena*120/log(2.0);
        }
        return(list(tp=throughput, featTP = feature_throughputs));
    }

    return(list(tp=throughput));
}

# Calculate throughput for a given pair of files.
TPpair <- function(filename1, filename2, fps = 120, pca = FALSE, amc = FALSE, res = FALSE, ftp = FALSE, war = FALSE, noise = 0) {
    a <- read.table(filename1);
    b <- read.table(filename2);

    # Eliminate features with zero variance in AMC format data, usually the
    # fingers of both hands.
    if (amc) {
        a$V34=NULL; a$V46=NULL;
        b$V34=NULL; b$V46=NULL;
    }

    datams <- remove_duplicate_frames(a, b);
    a <- datams[[1]];
    b <- datams[[2]];

    if (noise > 0) {
        a <- add_noise(a, noise);
        b <- add_noise(b, noise);
    }

    a <- normalize_features(a);
    b <- normalize_features(b);


    return(throughput(a, b, fps = fps, pca = pca, res = res, ftp = ftp, war = war));
}

# Calculate throughput for all aligned sequence pairs in the directory (by
# default the current working directory). Boolean parameters for the use of
# Principal Components analysis and AMC data (instead of coordinate data).
# Currently sequences 1 and 2, 3 and 4, 5 and 6 and so on are assumed to
# be repetitions of the same sequence.
#
# Framerate (fps) can be specified if it differs from 120.
#
# The function assumes that the original amc or coordinate files are in the
# working directory. (Note: all files should be for the same test subject!)
# Files should be named <sequence number>.txt, e.g. "15.txt".
#
# A subdirectory named 'aligneddata' should contain the aligned data
# in files that are named in the following way:
#         <seq1 #>_ali_<seq2 #>.txt
# Example: "14_ali_15.txt".
TPdir <- function(fps = 120, pca = FALSE, amc = FALSE, res = FALSE, ftp = FALSE, war = FALSE, noise = 0) {
    if (amc) {
        files <- dir(".", "^[[:digit:]]+.amc");
    } else
        files <- dir(".", "^[[:digit:]]+.txt");
    sequences <- length(files);
    
    total_throughputs <- array(0, c(sequences/2, 2));
    feature_throughputs <- array(list(NULL), c(sequences/2, 2));
	
    for (i in 1:(sequences/2)) {
        j = 2*i;
        k = j-1;
        file1 <- sprintf("aligneddata/%d_ali_%d.txt", k, j);
        file2 <- sprintf("aligneddata/%d_ali_%d.txt", j, k);
        total_throughputs[i,1] <- TPpair(file1, file2, fps = fps, pca = pca, amc = amc, res = res, ftp = ftp, war = war, noise = noise)$tp;
        total_throughputs[i,2] <- TPpair(file2, file1, fps = fps, pca = pca, amc = amc, res = res, ftp = ftp, war = war, noise = noise)$tp;
        if (ftp) {
            feature_throughputs[i,1] <- list(TPpair(file1, file2, fps = fps, pca = pca, amc = amc, res = res, ftp = ftp, war = war, noise = noise)$featTP);
            feature_throughputs[i,2] <- list(TPpair(file2, file1, fps = fps, pca = pca, amc = amc, res = res, ftp = ftp, war = war, noise = noise)$featTP);
        }
    }
    if (ftp) {
        return(list(total=total_throughputs, feature=feature_throughputs));
    } else {
        return(list(total=total_throughputs));
    }
}
