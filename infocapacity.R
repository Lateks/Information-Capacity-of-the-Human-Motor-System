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

twopie = 2*pi*exp(1);

# return stochastic complexity of sequence r with side information v
evalSC <- function(v, r) {
    n <- nrow(r);           # number of features
    totCL <- 0;             # total code-length over all features
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

        b <- qr.coef(Aplus, r[,i]);      # least squares fit
        res <- qr.resid(Aplus, r[, i]);  # residuals

        if (sum(res^2) < .01) print(c(i,sum(res^2)))  # warning about possibly too low variance

        rss <- sum(res^2)

        # Rissanen's classic two-part MDL code-length (feel free
        # to replace by more advanced universal code).
        totCL <- totCL + n/2*log(twopie*rss/n) + k/2*log(n);
    }
    return(totCL);
}

# return residuals from ar(2)
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
SI <- function(a,b) {
    n <- min(nrow(a), nrow(b))-2;
    totSIa <- 0;
    
    resa <- evalR(cbind(a[1:n,],a[1:n+1,]), a[1:n+2,]);
    resb <- evalR(cbind(b[1:n,],b[1:n+1,]), b[1:n+2,]);

    # Not sure if this is what we want to do.
    # This produces way too large results when summing
    # over all features.
    for (i in 1:length(resa)) {
        lm.ra = lm(resa[[i]] ~ resb[[i]]);
        resaf = fitted(lm.ra);

        rssa <- sum(resa[[i]]^2);
        rssaf <- sum(resaf^2);

        totSIa <- totSIa + (n/2*log(rssa/rssaf) - log(n)*3/2);
    }

    return(list(n=n, SIa=totSIa));
}

# complexity of a single (multivariate) sequence
SC <- function(a) {
    n <- nrow(a)-2;
    if (n < 1) return(list(n = 0,k = 1,cl = 0,rate = 0));
    cl <- evalSC(cbind(a[1:n,],a[1:n+1,]), a[1:n+2,]);
    return(list(n = n,cl = cl,rate = cl/n));
}

# complexity of a (multivariate) sequence given another
SCcond <- function(a,b) {
    n <- min(nrow(a), nrow(b))-2;
    cl <- evalSC(cbind(a[1:n,],a[1:n+1,],b[1:n+2,]), a[1:n+2,]);
    return(list(n=n,n,cl=cl,rate=cl/n));
}
   
# Process the given data matrices by removing rows duplicated in
# time warping and by normalizing features.
processdata <- function(a, b) {
    # remove rows duplicated in sequence 'a' due to alignment
    skip <- matrix(FALSE, nrow(a));
    skip <- rowSums((a[2:nrow(a),]-a[1:(nrow(a)-1),])^2) < 0.001;
    a <- a[!skip,];
    b <- b[!skip,];

    # normalize features
    a <- t(apply(a,1,'-',apply(a,2,mean)));
    a <- t(apply(a,1,'/',sqrt(apply(a,2,var))));
    a[is.nan(a)]<-0;
    b <- t(apply(b,1,'-',apply(b,2,mean)));
    b <- t(apply(b,1,'/',sqrt(apply(b,2,var))));
    b[is.nan(b)]<-0;

    return(list(a, b))
}

# Calculate throughput for a given pair of matrices.
throughput <- function(a, b, fps = 120, pca = FALSE, res = FALSE) {
    lena <- nrow(a);

    if (pca) {
        pcs <- prcomp(a, retx = TRUE, center = TRUE, scale. = TRUE);

        sum <- 0;
        evecs <- 0;
        threshold <- 0.9*sum(pcs$sdev**2);
        for (k in 1:length(pcs$sdev)) {
            sum <- sum + pcs$sdev[k]**2;
            evecs <- k;
            if (sum >= threshold)
                break;
        }
        eigenvec <- pcs$rotation[,1:evecs];
        a <- a %*% eigenvec;
        b <- b %*% eigenvec;
    }

    if (res) # Return shared information determined by residuals
        return(SI(a, b)$SIa/lena*120/log(2.0))

    # compute stochastic complexity with and without condition
    ca <- SC(a)$cl;
    cab <- SCcond(a, b)$cl;

    # calculate throughput
    tp <- (ca - cab)/lena*fps/log(2.0);

    return(tp);
}

# Calculate throughput for a given pair of files.
TPpair <- function(filename1, filename2, fps = 120, pca = FALSE, amc = FALSE, res = FALSE) {
    a <- read.table(filename1);
    b <- read.table(filename2);

    # Eliminate features with zero variance in AMC format data, usually the
    # fingers of both hands.
    if (amc) {
        a$V34=NULL; a$V46=NULL;
        b$V34=NULL; b$V46=NULL;
    }

    datams <- processdata(a, b);
    a <- datams[[1]]; b <- datams[[2]];

    return(throughput(a, b, fps = fps, pca = pca, res = res));
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
TPdir <- function(fps = 120, pca = FALSE, amc = FALSE, res = FALSE) {
    if (amc) { # count amc files
        files <- dir(".", "^[[:digit:]]+.amc");
    } else # count coordinate data files
        files <- dir(".", "^[[:digit:]]+.txt");
    seqs <- length(files);

    M <- matrix(0, seqs, seqs);

    for (i in seq(1, seqs-1, by=2)) {
        j = i+1;
        file1 <- sprintf("aligneddata/%d_ali_%d.txt", i, j);
        file2 <- sprintf("aligneddata/%d_ali_%d.txt", j, i);
        M[i,j] <- TPpair(file1, file2, fps = fps, pca = pca, amc = amc, res = res);
        M[j,i] <- TPpair(file2, file1, fps = fps, pca = pca, amc = amc, res = res);
    }

    return(M);
}
