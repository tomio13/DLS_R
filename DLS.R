#' Dynamic light scattering analysis functions, especially designed for
#' the ALV DLS system. However, the analysis and plotting tools are
#' useful for any other data sets.
#' They may require lists with specific elements, look into the
#' code and comments for further details.
#' Author: T. Haraszti (haraszti@dwi.rwth-aachen.de)
#' Date: 2018 -
#' Licence: CC-BY-4
#' Warranty: None


MSD <- function(a, N=3){
    #' Calculate he mean squared displacement assuming that the correlation function
    #' is exp( -MSD*q**2), thus estimate: -log(correlation)/q**2
    #' if q is in 1/micron, we get it in micron**2
    #' Kills NaN and infinite elements to -1
    #'
    #' @param a     a DLS object. A list containing 'correlation' and 'q' elements at least
    #'              The correlation element is a 2D array with first column the time, the others
    #'              are the normalized correlation functions exp(-MSD*q**2).
    #' @param N     Dimensionality, it should be 3 for DLS because we detect the diffusion in 3D
    #' @return  an array containing the time and the MSD values

    msd <- -log(as.array(a$correlation))*N/(a$q**2);
    #first column is tau, restore it!
    msd[,1] <- as.array(a$correlation[,1]);

    msd[is.infinite(msd) ] <- -1;
    msd[is.nan(msd)] <- -1;
    return(msd);
}


###############################################################################
# Helper functions for distributions

seq.log <- function(start, end, len=50){
    #' logaritmic sampled number sequence
    #' using seq(log(start), log(end), length=len)
    #'
    #' the resulted array is equidistant in log space
    #' @param start double  the start value
    #' @param end double    the end value
    #' @param len integer   length of the resulted array
    #'
    #' @return an 1D array including both ends

    return(exp(seq(log(start), log(end), length=len)))
}


smooth.dist <- function(dls.fit,
                        Rh.array= NULL,
                        norm=TRUE,
                        epsilon= 1E-3,
                        method='natural',
                        which= 'I',  ...){
    #' Calculate a smoothened version of the size distribution
    #' based on spline() but with logarithmic size array
    #' @param dls.fit   a fit list returned by analyse.size below
    #' @param Rh.array  the new size array where we want to have the
    #'                  points. Use seq.lod to generate one
    #' @param norm      should the distribution be normalized to
    #'                  its absolute maximum?
    #' @param epsilon   below this value, all interpolated dist values are set to zero
    #' @param method    passed to spline as method
    #' @param which     I,V or N for intensity, Volume or Number distribution
    #' @param ...       passed to spline for further parameters
    #'
    #' @return  a list containing x and y arrays

    if (is.null(Rh.array)){
        # Rh.array=exp(seq(log(0.01), log(1E4), len=250)),
        Rh.array=exp(seq(log(min(dls.fit$Rh.array)), log(max(dls.fit$Rh.array)), len=250))
    }

    x <- dls.fit$Rh.array;
    lx <- log(x);

    if (which == 'N' && !is.null(dls.fit$N.dist)){
        y <- dls.fit$N.dist
    } else if(which == 'V' && !is.null(dls.fit$V.dist)){
        y <- dls.fit$V.dist
    } else {
        y <- dls.fit$dist
    }
    x.out = log(Rh.array)

    ret <- spline(lx, y, method= method, ..., xout= x.out);
    clean.ret <- approx(lx, y, xout= x.out, yleft= 0, yright= 0, rule= 2);

    #ret <- predict(f.sp, log(Rh.array));
    #overwrite x with the original (not log) Rh.array
    ret$x <- Rh.array;
    ret$y[ret$y < 0] <- 0;
    ret$y[clean.ret$y < epsilon] <- 0;
    # data are normally normalized to sum(y) = 1
    # normalize to the maximum:
    if(norm == TRUE)   ret$y <- ret$y / max(abs(ret$y))

    return(ret);
}


################################################################################
# analysis functions: cumulant and size distribution


analysis.cumulant <- function(dls, channels=c(1), lin.fit.limit = -2, start=5){
    #' perform a cummulant analysis on the dls$correlation data
    #' this should be an array with first column the tau values in milliseconds,
    #' further columns are correlator channels (if there is more than one)
    #'
    #' @param dls       a list of dls data, e.g. from read.ALV
    #' @param channels  which channel to fit a single or an array of values (default 1)
    #' @param lin.fit.limit     the logarithm above this value will be
    #'                          taken (-1 --> 1/e) (default -2)
    #'
    #' @return
    #' a list containing the fitted curves and extracted information

    # kb is 10^-23 J/K; D is in micron^2/s thus 10^-12, mPas thus 10^-3: 23-15 = 8
    # 1 nm is 10^-9 m, thus 10^-8 = 10 nm... keep a 10x in kB, then we have Rh in nm
    kB <- 13.8064852; #J/K --> nm
    ncol <- length(channels);
    N <- nrow(dls$correlation);
    tau <- rep(dls$correlation[,1], ncol);
    y <- as.vector(dls$correlation[,-1][,channels]);
    # limit to start:N:
    tau <- tau[start:N];
    y <- y[start:N];

    if( is.null( dls$q ) ){
           q <- 4*pi*dls$n*sin(dls$theta/2*pi/180)/dls$lambda*1000; #scattering vector 1/micron
    } else {
        q <- dls$q;
    }
    # quick and dirty cumulant:
    if( min(y) < 0){amp <- 0;}
    else{           amp <- min(y);}

    # block log() from going crazy:
    yls <- y - amp;
    yls [ yls <= 0 ] <- exp(-10);
    yls <- log( yls );
    # limit the fit to the first part:
    indx.fit <- yls > lin.fit.limit;
    yls.fit <- yls[indx.fit];
    tau.fit <- tau[indx.fit];

    ord.indx <- order(tau.fit);
    yls.fit <- yls.fit[ord.indx];
    tau.fit <- tau.fit[ord.indx];

    # the fit is ln(g-1) = 2*(K0 - Gamma*tau + K2/2*tau^2-K3/(2*3)*tau^3)
    tau.fit.2 <- tau.fit^2/2; #in the fit we have K2*tau**2/2!
    tau.fit.3 <- -tau.fit^3/6; #in the fit we have K3*tau**3/3!
    tau.fit <- -tau.fit;

    # lin.fit <- lm( yls.fit ~ 1 + tau.fit + tau.fit.2);
    lin.fit <- lm( yls.fit ~ 1 + tau.fit + tau.fit.2 + tau.fit.3);

    # we are fitting g(1)**2, with 2x the coefficients of g(1)(tau).
    K.lin <- coef(lin.fit)[-1]/2;  #coeffs = 2*(K0,K1,K2,K3)
    # K2.lin <- coef(lin.fit)[3]/2;   #coeff3 = 2*K2
    # K3.lin <- coef(lin.fit)[4]/2;   #coeff4 = 2*K3
    names(K.lin) <- c('K1','K2','K3');
    D.lin.avg <- K.lin[1]/(q**2)*1000; #micron**2/second
    Rh.lin <- kB*dls$temperature/(6*pi*dls$eta*D.lin.avg); #in nm

    # now go for a nonlinear fitting:
    f <- function(p, x.array, y.array, weights=rep(1.0, length(x.array)) ){
        sum( weights*(p[1]*exp( -p[2]*x.array + p[3]*x.array^2 - p[4]*x.array^3) - y.array)^2 );
    }

    plot(-tau.fit, yls.fit, xlab= expression(paste(tau,', ms', sep='')),
         ylab=expression(paste('log(',g^(2),'(',tau,')-1)',sep=''))
    )
    lines(-tau.fit, predict(lin.fit), col='red', lwd=3);
    # print(lin.fit)

    # start <- c( exp(coef(lin.fit)[1]), K.lin[1], K.lin[2], K.lin[3] );
    start <- c( exp(coef(lin.fit)[1]), K.lin[1], 0, 0);
    cat('Start values:', start,'\n')

    fit <- optim( start, f,
                #method = 'L-BFGS-B',
                #lower = list(beta=0.1, k1=1E-8, k2= 0, k3= 0),
                #upper = list(beta= 1.2, k1 = 1, k2 = 0.1, k3 = 0.1),
                control = list(maxit = 5000), x.array = tau, y.array= y);
                #alternatively
                #control = list(maxit = 1000), x.array = tau, y.array= y, weight);

    # here we did not modify the coefficients
    parms <- fit$par
    names(parms) <- c('beta','K1','K2','K3');

    K1 <- parms['K1']/2; #coef1 = 2 K1
    K2 <- parms['K2'];  #coef2 = 2 K2/2
    K3 <- 3*parms['K3'];    #coef3 = 2* K3/6

    D <- K1 / q^2 * 1000;
    names(D) <- "diffusion coeff";
    Rh <- kB*dls$temperature/(6*pi*dls$eta*D); #in nm
    names(Rh) <- "Hydrodyn Radius";

    # for return purposes:
    cf <- coef(lin.fit);
    y.lin.fitted <- exp( cf[1] - tau*cf[2] + cf[3]/2*tau^2 - cf[4]/3*tau^3) + amp;
    y.fitted <- parms[1]*exp(-parms[2]*tau + parms[3]*tau^2 - parms[4]*tau^3);

    lines( tau, log(y.fitted), col='green', lwd=3)

    return(list(
                lin.fit = lin.fit, q=q,
                #P = K2/K1^2
                P.lin = K.lin[2]/K.lin[1]^2,
                P = K2/K1^2,
                lin.fit.indx=indx.fit,
                K.lin = K.lin,
                K = c(K1, K2, K3),
                fit = fit,
                par = parms,
                D = D, Rh = Rh,
                D.lin.avg= D.lin.avg, Rh.lin= Rh.lin, lin.min= amp,
                y.lin.fitted = y.lin.fitted,
                y.fitted = y.fitted,
                tau = tau
                ));
}


analysis.sizes <- function(dls, Rh.array= seq.log(8.66264E-3, 8.66264E3, len=50),
                           g2.channel = c(2),
                           tau.min= -1, tau.max= -1,
                           start.parm = NULL,
                           # weight is typically in the order of 0.3 ... 0.6, so do not top far
                           # above it:
                           weight.max = 0.66,
                           #alpha.array = exp( seq( log(5E-3), log(50), len=30))){
                           alpha.array = seq.log(50, 5E-3, len=50),
                           factor= 0.99,
                           inherit = TRUE,
                           plot.fit = TRUE){
#'  A regularized fitting to construct the correlation function in dls as a sum of exponentials
#'  The method is between the CONTIN and a direct regularized fit
#'  It fits the g(2)-1 function (squared) not its square root as in the literature, avoiding
#'  the extra squaring of the result.
#'  Using boxed constrains, the parameters are limited to the [0,1] closed interval.
#'  The best alpha is searched in the alpha.array, to provide the minimum of the fit error f.norm
#'  If inherit is TRUE, the previous result will be fed as starting parameters
#'  If start.parm is NULL, a zero array is used for starting parameter, with one point set to 0.5
#'  for the radius corresponding to where the correlation function falls to 1/e it amplitude
#'
#'  Running with decreasing alpha values, we get the most details of the fit which still
#'  improved the fit by factor to the previous best.
#'  (High alpha means smoothness dominates the distribution.)
#'
#'  The algorithm was inspired by:
#"      A. Scotti et al. 'The CONTING algorithm and its application to determine
#'      the size distribution of microgel suspensions' in The Journal of Chemical Physics
#'      vol. 142, 234905 (2015)
#' and the original CONTIN papers:
#'      S. W. Provencher, 'CONTIN: A general purpose constrained regularization program
#'      for inverting noisy linear algebric and integral equations' in
#'      Computer Physics Communications vol. 27, 229 - 242 (1982)
#' and
#'      S. W. Provencher, 'A constrained regularization method for inverting data represented
#'      by linear algebric or integral equations' in
#'      Computer Physics Communications vol.: 27, 213 - 227 (1982)
#'
#'
#' @param dls           the dls list containing the correlation data as well as
#'                      parameters, such as scattering vector, angle, temperature...
#' @param Rh.array      an array of radii in nanometers, nm
#'                      Select the radius distribution carefully. Too low walues will bias
#'                      the plateau. The large range may have to go to up to non-physical
#'                      ranges, e.g. more than 50 microns.
#'                      Default array is taken from the ALV software
#'
#' @param g2.channel    name of the channel in the correlation array
#' @param tau.min       start fitting from this delay value (ms)
#' @param tau.max       fit up to this delay time (ms)
#'                      if they are set to -1, then the [3:length(g2)] indices are taken
#'                      (ALV standard values)
#' @param start.param   the distribution array to start with. By default it is estimated
#'                      using the exponential decay of the g1 function
#' @param weight.max    maximal weight allowed. It may influence the convergence, because
#'                      g2 values close to zero may kick it up errorneously
#'
#' @param alpha.array   an array of alpha values to be scanned through. It is good to have
#'                      it logarithmically spaced
#' @param factor        a new error is suitable new result if it is less than factor*old error
#'                      if the alphas are increasing then less or equal is taken
#' @param inherit       Feed the best fit to the start values of the next alpha fit
#' @param plot.fit      a Boolean, if set, plot the correlation, fit and residuals in a split plot

    if(length(Rh.array) < 1){
        cat("Invalid size array Rh\n");
        return;
    }
    # change display precision:
    old.opt <- options();
    options(digits = 4, scipen=-2);

    kB <- 13.8064852; #J/K --> nm
    N <- length(Rh.array);
    D.array <- kB*dls$temperature / (6*pi*dls$eta*Rh.array); #in micron**2/second
    # from:  D.lin.avg <- K.lin[1]/(q**2)*1000; #micron**2/second; already divided by 2
    # remove the 2* factor to work on the g(1)(tau):

    # unit matching to G in 1/ms = kHz
    # G.array <- dls$q^2*D.array/1000;
    G.array <- 2*dls$q^2*D.array/1000;

    # tau <- dls$correlation[,'tau'];
    tau <- dls$correlation[,'tau'];
    tau.range = range(tau) # to be used in plotting

    if(tau.min <= 0){
        start <- 3;
    } else {
        start <- max(which(tau < tau.min));
        if(start < 1)          start <- 1;
    }

    if(tau.max <= 0){
        # N.tau <- 129; #default at AVL
        N.tau <- length(tau)
    } else {
        N.tau <- max(which(tau < tau.max))
        if(N.tau < start )     N.tau <- length(tau);
    }
    cat('Fit starts at:',start, tau[start],'ms \n');
    cat('Fit ends at:',N.tau, tau[N.tau],'ms \n');
    tau <- tau[start:N.tau];

    func.array <- exp(outer(tau, -G.array, "*"));
    # now each column is the time, the rows are the Rh.array related decay rate (Gamma)
    # reduce the matrix to something useful
    # and our result transforms from g2 to y:
    g2 <- dls$correlation[start:N.tau,g2.channel]
    # backup the original
    g2.orig <- g2
    # sometimes the function goes off, it should be up to
    # 1 (or 1+max(g2) for other correlators)

    # readjust because start > 0 and end < N
    N.tau <- length(g2);
    # background should be about 0:
    g2.bg <- mean(g2[(N.tau-3):N.tau])

    if(g2.bg < 0){
        g2.bg <- 0;
    }

    g2 <- g2 - g2.bg;

    if (weight.max > 0) {
        # weight from the maximized entropy paper:
        #   S-L Nyeao and B. Chu, Macromolecules 22:3998-4009 (1989)
        # formula (20) using B=1
        # weight <- (1 + g2)/(4 * g2);
        weight <- 1/(4 * g2) + 0.25
        # cancel oversized weights and those in the noise range:
        # (it is clear that this should not happen...)
        weight[weight < 1E-3] <- 0;
        # weigh.max used to be 500
        weight[weight > weight.max] <- weight.max
        # debug
        # print(weight)
    } else {
        # turn off weights
        weight <- rep(1, length(g2))
    }

    # normalize to c.a. 1:
    # g2.m <- max(g2)
    # use the plateu (hopefully) of the first few points:
    g2.m <-  mean(g2[1:5])
    cat('normalizing to:', g2.m, '\n')

    # g2 <- sqrt(g2/g2.m)
    g2 <- g2/g2.m
    g2[is.nan(g2)] <- 0

    # the CONTIN algorithm fits the correlation function with
    # the linear combination of a set of sizes.
    # However, since this is a mathematically ill defined situation,
    # it employs a penalty on how compact the resulted size distribution
    # is, using the second derivative of this distribution.
    #
    # create a second derivative operator matrix:
    Om <- matrix(0, ncol=N, nrow=N);
    diag(Om) <- -2;
    indx.1 <- matrix(1:N, ncol=2, nrow=N);
    indx.2 <- indx.1;
    indx.1[,1] <- indx.1[,1] - 1;
    indx.2[,2] <- indx.2[,2] - 1;
    Om[indx.1[-1,] ] <- 1;
    Om[indx.2[-1,] ] <- 1;
    # multiplying this with the size distribution in x will give
    # the second derivative function. The last points have a tailing
    # error, which we drop.

    # solve the problem, we have then the error defined as:
    f <- function(x, alpha.parm=0.0, weights=rep(1, length(g2))){
            #delete the last and first element of the derivative, because it biases the result!
            sum((weights*(func.array%*%x - g2))^2) + sum((alpha.parm*Om%*%x)[-c(1,N)]^2)
    }

    # make an iterative solution to find the best stabilizer value
    # alpha should be an array
    if(length(alpha.array) > 0){
        alphas <- alpha.array
    } else {
        cat("Invalid alpha array, falling back to default\n");
        # decreasing means we start with something smooth, and allow for more and
        # more details if they improve the fit
        alphas <- seq.log(20, 5E-3, len=30);
    }

    if(alphas[1] > alphas[length(alphas)]){
        alphas.decreasing <- TRUE
    }
    else{
        alphas.decreasing <- FALSE
    }

    errs <- rep(0.0, length(alphas));
    x.norm <- rep(0.0, length(alphas));
    f.norm <- rep(0.0, length(alphas));

    dist.array <- array(0, dim= c(length(Rh.array), length(alphas)))

    i <- 1;
    # the error will never be so high...
    err.min <- Inf;
    end.alpha <- -1;
    end.fit <- NULL;


    cat('running on', length(Rh.array),'size values\n');
    cat('and',length(tau),'data points\n');
    # start with some reasonable parameters
    # assuming a single exponential curve, one can estimate the decay rate
    #  at the 1/e time value, where tau/t = 1
    # this is the Gamma*tau = 1, then Gamma is related to an R
    # this single peak has some value, e.g. 0.5 to start with
    indx.tau.R <- max(which(g2 > exp(-1)))
    tau.R <- tau[indx.tau.R]
    R <- dls$q**2*kB*dls$temperature*tau.R / (3000*pi*dls$eta);
    cat('1/e time is:', tau.R,'value:', g2[indx.tau.R], '\n');
    cat('Estimated radius from exponential decrease is', R, 'nm\n');

    set.parm <- TRUE;
    if(! is.null(start.parm)){
        if(length(start.parm) == length(Rh.array)){
            par.start <- start.parm;
            set.parm <- FALSE;
        }else{
            cat('Start parameters have to have a length of', length(Rh.array),'\n');
            set.parm <- TRUE;
        }
    }
    if(set.parm == TRUE){
        cat('Setting up start paramters for fitting automatically\n');
        par.start <- rep(0, N);
        #cat('Setting R:',R,'\n');
        cat('Setting R:',R,'nm to 0.5\n');
        par.start[max(which(Rh.array <= R)) ] <- 0.5
    }

    # run through the reularizer series:
    for(alpha in alphas){
        # these start values give a quite good start estimate:
        #  the end values are the amplitudes the various exponents
        #  building up the sum forming g2(tau)
        # Constrains are now built into the error landscape...
        fit <- optim(par.start, f, method='L-BFGS-B',
                     upper= rep(1.0, length(par.start)),
                     lower= rep(0.0, length(par.start)),
                     control=list(maxit= 500),
                    alpha.parm = alpha, weights= weight);

        # cat('i:',i,'alpha:', alpha, "error", fit$value,'\n', sep='\t');
        errs[i] <- fit$value; #the value of f(), which is the estimated error
        x.norm[i] <-  sum((Om%*%fit$par)[-c(1,N)]^2);

        # the error is recalculated without weight and alpha, only the
        # exponential distribution:
        f.norm[i] <- f(fit$par, alpha.parm=0, weights= weight);
        dist.array[,i] <- fit$par;

        # while it is quite bad for the constrOptim, it works with optim
        # update the start with the result of the previous run
        # it causes roughening only when needed (or so we hope)
        if(inherit == TRUE)        par.start <- fit$par;

        # instead of running it again, check if it is the best
        # so far and store it if it was so
        cat('i:',i,'alpha:', alpha, '    error', fit$value,'    chi2',f.norm[i],'\n', sep='\t');
        # while for small alpha it gives the same result, it should be sane
        # following the chi2 minima in the game: which is then the best fit?
        # originally we used the complete error to estimate improvements
        # we can use the fit error for this purpose
        # we use < here so the smoothest result changes only if the
        # new features improve the fit!
        #   Here the only difference is that if alphas increase, we allow '=' too
        #   this forces the smoother from two equivalent fits for that case
        if(f.norm[i] < factor*err.min ||
           (alphas.decreasing == FALSE && f.norm[i] <= factor*err.min)){
            #err.min <- fit$value;
            err.min <- f.norm[i];
            end.fit <- fit;
            end.alpha <- alpha;
        }
        #next:
        i<- i +1;
    }

    # transform the fitted data to the original g2
    # put back to the square, for the similarity to the data.
    # y.fitted <- g2.m*(func.array%*%end.fit$par)**2 +g2.bg;
    y.fitted <- g2.m*(func.array%*%end.fit$par) +g2.bg;

    # print some summary:
    cat('Best alpha was:', end.alpha,'\n', sep='\t');
    cat('Best error:', err.min,'\n', sep='\t');
#    cat('Summary of errors\n');
#    err.matrix <- matrix( c(alphas, errs, x.norm, x.norm*alphas, f.norm), ncol=5)
#    colnames(err.matrix) <- c('alphas','error','reg.error','reg*alpha','funciton error');
#    print(err.matrix)

    # plot.DLS(dls, ch= 'ch1');
    # lines(tau, y.fitted, col='red', lwd=2);
    # indicate the noise limit defined by tau.max:
    # abline( h = g2[N.tau], col='blue', lwd=2);

    # How to get other forms of distribution?
    # form factor: p.R
    # theoretically dist = N.dist * I.R = N.dist * V**2 * p.R
    # Rh.array is in nm, q is in 1/micron:
    qR <- Rh.array/1000.0 * dls$q

    # isotropic form factor of a sphere (see: http://gisaxs.com/index.php/Form_Factor:Sphere)
    # for a sphere:
    # I.R = [3 \Delta\rho V (sin(qR) - qR cos(qR))/(qR)^3]^2, and V = 4 pi R^3 / 3
    # from this for isotropic form factor intensity one gets:
    # I.R = \Delta\rho^2 (4 \pi)^3 ((sin(qR) - qR cos(qR))/q^3)^2
    #
    # To get a number distribution, we need to multiply the intensity distribution with
    # I.R
    # however, we have no idea about the contrast factor (in polarisability)
    # thus we have only a proportionality part here...
    # Because of this, we can actually drop the prescaler (4 \pi)^3
    # we can renormalize the distribution to sum = 1 instead if it was needed...
    I.R <- ((sin(qR) - qR*cos(qR))/ dls$q**3)**2

    dist <- end.fit$par
    N.dist <- dist/I.R

    # V.dist = N.dist * V = dist/(V p.R)
    V.dist <- N.dist * 4/3*pi*(Rh.array/1000.0)**3

    # renormalize distributions using sums,
    # because integral would be a bit tricky in logarithmic
    # R space....
    N.dist <- N.dist/sum(N.dist)
    V.dist <- V.dist/sum(V.dist)

    # recover display precision
    options(old.opt);

    if (plot.fit) {
        # prepare a plot
        p <- par(c('mar', 'fig', 'mfrow'))
        par(mfrow= c(2,1), fig=c(0, 1, 0.3, 1), mar= c(0, 4.5, 2, 1), xaxt='n')
        # plot the original correlation curve and the fit:
        plot.DLS(dls, ch= 'ch1', dual= FALSE)
        lines(tau, y.fitted, col='blue', lwd=2)
        abline(h= 0, col='green')

        # now, the residues;
        par(new= TRUE, fig= c(0, 1, 0.0, 0.3), mar= c(5, 4.5, 0, 1), xaxt='s')
        # we restore the g2 to the original channel, but without the cut part
        plot(tau, g2.orig - y.fitted,
             log='x',
             type='l',
             lwd= 2,
             # ylab= expression(paste('log(',g^(2),'(',tau,')-1) error',sep='')),
             ylab= 'error',
             xlab= expression(paste(tau,', ms', sep='')),
             xlim = tau.range
        )
        abline(h= 0, col='green')
        # restore parameters
        par(p)
    }

    return(list(fit= end.fit,
                # also export the time stamp
                filename = dls$filename,
                time = dls$time,
                temperature= dls$temperature,
                eta = dls$eta,
                q = dls$q,
                #fit.0 = fit.0,
                tau = tau,
                fitted = y.fitted, weight= weight, residuals= g2.orig - y.fitted,
                g2= g2.orig, g2.norm = g2.m, g2.bg = g2.bg,
                D.array = D.array,
                G.array = G.array,
                Rh.array = Rh.array,
                dist = dist,
                dist.norm= dist / max(dist),
                N.dist = N.dist,
                V.dist = V.dist,
                x.norm = x.norm, f.norm = f.norm,
                alphas = alphas, end.alpha= end.alpha, errs = errs,
                dist.array = dist.array));
}


################################################################################
# plot functions designed for DLS


plot.correlation <- function(a, log='x', add= FALSE, ...){
    #' plot the correlation funciton from a DLS measurement
    #' expects a 2 or 4 channel data set
    #' Parameters:
    #' an array, where first column is tau
    #' further columns are the correlation function
    #' return value: none
    tau <- a[,1]
    if( dim(a)[2] < 2 ){
        cat(' array should have more than one column\n');
        return
    }

    if (add == TRUE) {
        points(a[,1], a[,2], ...)
    } else {
        plot(a[,1], a[,2], xlab='delay time, ms', ylab='correlation', log= log, ...)
    }

    M <- dim(a)[2]
    if( M > 2){
        for( i in  3:M) {
            if( sum(a[,i]) > 0) points(a[,1], a[,i], col=i-1)
        }
    }
}


plot.DLS <- function(dls, fit=NULL, ch=c('ch1', 'ch2'), dual= TRUE, log='x', xlim= NULL,...){
    #' Quick plotting of the autocorrelation in semilog format as it is common
    #' if fit is provided for the whole tau range, add the prediction of the fit
    #' @param dls   a DLS imported data list
    #' @param fit a cumulative fit list, not a fit object! (no predict() method!)
    #' @param ch character array    which channels to plot
    #' @param dual Boolean  if true, plot the residuals on the bottom
    #' @param log character string   which axis should be a log plot
    #' @param ... passed to plot()

    if (!is.null(fit) && dual){
        # we want a double plot
        p <- par(c('mar', 'fig', 'mfrow'))
        par(mfrow= c(2,1), fig=c(0, 1, 0.3, 1), mar= c(0, 4.5, 2, 1), xaxt='n')
    }
    firstplot <- TRUE
    tau <- dls$correlation[,'tau']
    for (this.ch in ch) {
        if (firstplot) {
            plot(tau, dls$correlation[,this.ch],
                 log= log,
                 xlab = "delay time, ms",
                ylab= expression(paste(plain(g)^2,'(',tau,')', sep='')),
                xlim= xlim,
                ...
            )
            firstplot <- FALSE
        }
        else {
            points(tau, dls$correlation[, this.ch], col='blue')
        }
    }

    if(!is.null(fit)){
        prediction <- fit$par[1]*exp(-fit$par[2]*tau + fit$par[3]*tau^2 -fit$par[4]*tau^3 )
        # complete the top plot
        lines(tau, prediction,
            col='red', lwd=2);

        if(dual) {
            # and make the bottom plot
            par(new= TRUE, fig= c(0, 1, 0.0, 0.3), mar= c(5, 4.5, 0, 1), xaxt='s')
            plot(tau, dls$correlation[, 1] - prediction,
                type= 'l',
                xlab= 'delay time, ms',
                ylab= 'error',
                log = log,
                xlim = xlim,
               ...
            )
            abline(h= 0, col='green')
            par <- p
        }
    }
}


plot.DLS.dist <- function(ft,
                          which= 'I',
                          smooth= TRUE,
                          add= FALSE,
                          xlab="hydrodynamic radius, nm",
                          ylab= "intensity weight",
                          type='p',
                          log= 'x',
                          main ='',
                          col=1,
                          xlim= NULL,
                          ylim= NULL,
                          ...){
    #' shortcut to plot the size distribution from a fit resulted in analysis.sizes
    #' @param ft a fit list from a size analysis
    #' @param which     a character I, N or V for the type of distribution
    #' @param smooth    bool,   also plot a smoothened line across
    #'                          use smooth.dist above
    #' @param add       bool,   add to the plot, not a new plot
    #' @param xlab, ylab, type, log, main, col: plot parameters
    #'  all further parameters go to smooth.dist
    #'
    #' @return invisible return the smoothed distribution

    if (which == 'V' && !is.null(ft$V.dist)){
        y <- ft$V.dist
        ylab <- 'volume distribution'
    } else if (which == 'N' && !is.null(ft$N.dist)){
        y <- ft$N.dist
        ylab <- 'number distribution'
    } else {
        y <- ft$dist
    }

    y <- y/sum(y)
    #norm <- sum(ft$fit$par)
    #y <- ft$fit$par / norm

    if (add) {
        points(ft$Rh.array, y, col=col)
    } else {
        plot(ft$Rh.array, y, log= log,
            type= type,
            xlab= xlab,
            ylab= ylab,
            main=main,
            xlim= xlim,
            ylim= ylim,
            col=col)
    }

    sd = smooth.dist(ft, which=which, ...)
    if(smooth == TRUE){
        lines(sd$x, sd$y*max(y), col= col)
    }
    invisible(sd)
}


###############################################################################
# plot functions for various series
# these assume lists generated either the way the read.ALV does,
# or the analysis.sizes does.


plot.corr.series <- function(namelst, norm= TRUE, norm.range=c(4:8), legends= NULL, ...) {
    #' take a list of variable names,
    #' and plot the correlation function of each
    #' making time the legend
    #' For example having experiments loaded as names
    #' var.001, var.002, var.003...
    #' you can use lst(pattern='^var\\.') for namelst.
    #'
    #' @param namelst a list of names
    #' @param norm  bool    if true, normalize the correlation functions to 1
    #' @param norm.range array  which points to use for normalization
    #'
    #' @return nothing
    N <- length(namelst)

    for (i in 1:N){
        a <- get(namelst[i], env= .GlobalEnv)
        x <- a$correlation[,1]
        y <- a$correlation[,2]
        if (norm) {
            y <- y / mean(y[norm.range])
        }
        if(i == 1) {
            plot(x, y, log='x',
                 xlab= expression(paste(tau,', ms', sep='')),
                 ylab= expression(paste('g'^2, '(', tau, ')', sep='')),
                 col= 1, pch= i,
                 type = 'o', ...)
        } else {
            points(x, y, type='o', pch= i, col=i)
        }
    }

    if (!is.null(legends)){
        legend('topright', col=1:N, pch= 1:N, legend= legends)
    }
}


plot.dist.series <- function(namelst, Rh.array.refined= NULL, legends= NULL, ...) {
    #' take a list of variable names for fits made using analysis.sizes,
    #' and plot the size distribution from them
    #'
    #' @param namelst   an array of variable names containing fits
    #' @param Rh.arrah.refined  an array of radii to be used (optional)
    #' @param legend    array of strings to be used as legends
    #' @param ...       extra plot parameters
    #'
    #' @return nothing

    N <- length(namelst)
    for ( i in 1:N) {
        a <- get(namelst[[i]], env= .GlobalEnv)
        if (i==1) {
             plot(a$Rh.array, a$dist, log='x',
                 xlab='hydrodynamic radius, nm', ylab='intensity weight',
                 col= i,
                 pch= i, ...)
        } else {
            points(a$Rh.array, a$dist, col=i, pch= i)
        }
        sd = smooth.dist(a, Rh.array= Rh.array.refined, norm=FALSE)
        lines(sd, col=i, lwd=2)
    }

    if (!is.null(legends)) {
        legend('topright', col=1:N, pch= 1:N, legend= legends)
    }
}

###############################################################################
# functions to manipulate size distribution / correlation


g.2 <- function(tau, G.array, dist) {
    #' calculate the correlation function for DLS
    #' based on the delay times and an array for
    #' the delay times with their weights
    #' ideally the weights have a sum of 1
    #' @param tau       array of delay times (ms)
    #' @param G.array   array of exponent decay times (1/ms)
    #' @param dist      weigth of the decay times
    #'
    #' @return an array for every tau value

    func.array <- exp(outer(tau, -G.array))
    return(func.array %*% dist)
}


correlation.limit.sizes <- function(fit, Rh.min= NULL, Rh.max= NULL) {
    #' take a size fit object (list resulted by analysis.sizes)
    #' plot the size distribution and activate a selector
    #' so the user can select a range to be kept
    #'
    #' If Rh.min and / or Rh.max are specified, skip the selection plot,
    #' use them to define the range
    #'
    #' Calculate the correlation function of the selection,
    #'  and subtract from the original correlation function
    #' plot what is remaining and return it as a table
    #' Required fields in fit:
    #' g2 the correlation function
    #' tau the delay times
    #' g2.norm --> normalization factor for g2
    #' g2.bg --> background constant in the g2 data
    #' during fit the (g2 - g2.bg)/g2.norm was fit
    #'
    #' @param fit   a fit list object from analysis.sizes
    #' @param Rh.min    beginning of R-range to be kept, default NULL
    #' @param Rh.max    end of R-range to be kept, default NULL
    #'
    #' @return  an array containing delay times and correlation values

    if(is.null(fit$g2) || is.null(fit$tau) || is.null(fit$Rh.array) ||
       is.null(fit$dist) || is.null(fit$temperature) || is.null(fit$eta) ||
       is.null(fit$q)) {
        cat('some critical data are missing\n')
        cat('function requires fit having arrays for: g2, tau, Rh.array and dist at least\n')
        return()
    }

    if (is.null(Rh.min) && is.null(Rh.max)) {
        plot(fit$Rh.array, fit$dist, log='x', type='o',
             main='click on the start and end area of distribution to be kept',
             sub='right click when done'
             )
        selected <- locator()
        # print(selected)
        Rh.range <- selected$x
        if (length(Rh.range) >= 2) {
            Rh.range <- tail(selected$x, 2)
        } else {
            cat('one or no points were selected\n')
            return()
        }
        Rh.min <- min(Rh.range)
        Rh.max <- max(Rh.range)
        # cat('selected:', Rh.min, Rh.max, '\n')
    }

    if (is.null(Rh.min)) {
        Rh.min <- min(fit$Rh.array)
    }
    if (is.null(Rh.max)) {
        Rh.max <- max(fit$Rh.array)
    }

    if (Rh.min == Rh.max) {
        cat('empty range is selected\n')
        return()
    }

    indx <- (fit$Rh.array >= Rh.min) & (fit$Rh.array <= Rh.max)
    Rh.array <- fit$Rh.array[indx]
    dist <- fit$dist[indx]
    points(Rh.array, dist, col='blue', pch= 16, cex= 1.1)

    kB <- 13.8064852; #J/K --> nm
    g2.norm <- ifelse(is.null(fit$g2.norm), 1, fit$g2.norm)
    g2.bg <- ifelse(is.null(fit$g2.bg), 0, fit$g2.bg)
    D.array <- kB*fit$temperature / (6*pi*fit$eta*Rh.array); # in micron^2 / sec.
    # unit matching to G in 1/ms = kHz
    G.array <- 2*fit$q^2*D.array/1000;
    func.array <- exp(outer(fit$tau, -G.array))
    # calculate the correlation but convert it to the experimental curve
    # scaling it and applying a constand background shift based on the known fit
    g2.calc <- g2.norm * func.array %*% dist + g2.bg
    g2.diff <- fit$g2 - g2.calc

    dev.new()
    plot(fit$tau, g2.diff, log='x', type='o',
         xlab=expression(paste(tau, ', ms', sep='')),
         ylab= expression(paste('g'^'(2)', '(', tau, ')', sep='')),
         main='residual g2'
         )

    correlation <- matrix(c(fit$tau, g2.diff), ncol= 2)
    colnames(correlation) <- c('tau', 'ch1')

    return(list(
                      g2.orig = fit$g2,
                      g2.calc = g2.calc,
                      correlation= correlation,
                      q = fit$q,
                      eta = fit$eta,
                      temperature = fit$temperature,
                      D.array = D.array,
                      G.array = G.array,
                      Rh.array = Rh.array,
                      dist = dist
             )
    )
}
