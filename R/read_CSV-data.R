read.DLS.zetasizer <- function(filename, theta=90.0,
                         temperature= 25.0,
                         eta = 0.89,
                         n = 1.334,
                         lambda= 632.8,
                         microsec= TRUE,
                         sep= ',',
                         dec= '.',
                         ch.names = c('ch1','ch2','ch3')
                         ) {
    #' read a csv file containing the correlation data from Zetasizer
    #' @details
    #' Read a csv file containing the correlation data from Zetasizer
    #' and form a list for the DLS processing, adding the missing parameters.
    #' The array is typically three X,Y columns containing the correlation
    #' delay time as x in microseconds, and g2-1 for correlation functions,
    #' for the steady state, transient and unfiltered versions.
    #' However, this may be somewhat different how the export was performed.
    #'
    #' @param filename  text, the file to read
    #' @param theta     float, the scattering angle in degrees
    #' @param temperature float, in degrees Celsius
    #' @param eta       float, the viscosity in mPas (cP)
    #' @param n         float, the refractive index of the medium
    #' @param lambda    float, the wavelength of the laser in nm
    #' @param microsec  Boolean, if the time is in microseconds
    #'                  instead of milliseconds
    #' @param sep       character, the field separator
    #' @param dec       character, decimal separator, i.e. '.' or ','
    #' @param ch.names a list of 3 names for the correlation names
    #'
    #' @return      a list containing all the information
    #'              time, filename, lambda, temperature,
    #'              n, eta, theta, theta.rad, q, correlation
    #' @export

    t.names <- c('tau', 'tau.trans', 'tau.unfilt')

    a <- read.table(filename, as.is= TRUE, sep= sep, dec= dec,
                skip= 1, col.names = array(t(array(c(t.names, ch.names), dim=c(3,2))), dim=c(6)) )

    print(str(a))

    if (microsec == TRUE) {
        a[,c(1,3,5)] <- a[,c(1,3,5)]/1000.0
    }

    if (all(a[,1] == a[,3])) {
            a <- a[,c('tau', ch.names)]
    }
    result <- list()
    result$time <- 0;
    result$filename = filename;
    result$lambda <- lambda #in nm
    result$temperature <- temperature + 273.15 #in Kelvin
    result$n <- n #water
    result$eta <- eta #cP, water
    result$theta <- theta
    result$theta.rad <- pi*result$theta/180. #theta in radian
    result$q <- 4*pi*result$n*sin(result$theta.rad/2)*1000.0/result$lambda #1/micron
    result$correlation <- a
    return(result)
}
