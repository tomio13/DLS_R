read.ALV <- function(filename){
    #' read ALV correlation files
    #' @details
    #' Read the ASC file output of the ALV DLS system, extract
    #' metadata related to the experiment and the autocorrelation
    #' tables.
    #' Important metadata include date, time, sample name, temperature,
    #' medium refractive index, viscosity, scattering angle
    #'
    #' The file has a relatively fixed structure, which is read out in as text
    #' then processed according to a predefined schema (a set of fields expected)
    #'
    #' @param filename text, name of the file to read (with path)
    #;
    #' @return a list containing all information, also adding parameters such as
    #'          q (scattering vector), theta.rad the scattering angle in radians
    #'
    #' @export

    if( ! file.exists(filename) ){
        cat("Invalid filename ", filename, "is provided\n");
        return()
    }
    # we use a predefined set of fields which come line-by-line
    # in the file header
    fields <- c('date', 'time','sample',
                paste('sample', 1:10, sep='.'),
                'temperature', 'eta', 'n', 'lambda',
                'theta', 'duration', 'duration.float', 'time.stop', 'runs',
                'mode', 'mean.cr.0', 'mean.cr.1', 'mean.cr.2', 'mean.cr.3');

    cat('reading file', filename,'\n');
    a <- readLines(filename, skipNul=T, encoding='latin-1');
    empty_indx <- which(a == '');

    result <- list();
    result$filename <- filename;
    result$version <- a[1];

    i <- 2;
    for(t in fields){
        # cat('analyzing:', a[i], '\n')

        # these lines have a : and some \t in them
        # strip off \" signs and tabs, extra spaced up
        # and after the ':'
        # set useBytes = TRUE to avoid encoding trouble we may face
        # the device runs Win 7 with whatever MS encoding...
        # others try using UTF-8
        s <- gsub("(^.+:\\s+)",'',a[i], perl=T, useBytes=TRUE);
        # remove quotes
        s <- gsub("\"","",s);

        # cat('expected field:', t, '\n')
        # cat('processed line:', s, '\n')

        # for these we expect numbers, like temperature, refractive index, etc.
        if(i > 14 && i != 24){
            result[[t]] <- as.numeric(s);
        } else {
            result[[t]] <- s;
        }
        i <- i+1;
    }
    # interpret date time string for R:
    result$date <- strptime(result$date, format='%d.%m.%Y');
    result$time <- strptime(paste(result$date, result$time),
                            format='%Y-%m-%d %H:%M:%S');

    # merge the sample description lines to one string
    indx <- grepl('sample', names(result))
    s <- paste(
               trimws(
                      result[grepl('sample',names(result))]
                      ),
               collapse=' ')

    cat('sample:', s, '\n')
    result[indx] <- NULL
    result$sample <- s

    i <- i+1;

    # the next block of lines until an empty one should be
    # the correlation table
    if(a[i] == '"Correlation"'){
        i <- i+1;
        i1 <- empty_indx[(empty_indx > i)][1];
        result$correlation <- suppressWarnings(
                                              as.numeric(
                                                         unlist(
                                                                strsplit(a[i:i1], '\t')
                                                                )
                                                         )
                                              );
        result$correlation <- matrix(result$correlation, byrow=T, ncol=5);
        i <- i1+1;
        colnames(result$correlation)<-c('tau','ch1','ch2','ch3','ch4');
    }

    if(a[i] == '"Count Rate"'){
        i <- i+1;
        i1 <- empty_indx[(empty_indx > i)][1];
        result$count.rate <- suppressWarnings(
                                             as.numeric(
                                                        unlist(
                                                               strsplit(a[i:i1], '\t')
                                                               )
                                                        )
                                             );
        result$count.rate <- matrix(result$count.rate, byrow=T, ncol=5);
        i <- i1+1;
        colnames(result$count.rate)<-c('tau','ch1','ch2','ch3','ch4');
    }

    # diode count is proportional to laser power
    result$laser.power <- suppressWarnings(
                                          as.numeric(
                                                     unlist(strsplit(a[i], '\t')
                                                            )[2]
                                                     )
                                          );
    i <- i+1;

    result$theta.rad <- result$theta*pi/180;
    # scattering vector 1/micron
    result$q <- 4*pi*result$n*sin(result$theta.rad/2)/result$lambda*1000;
    # next is the cumulant results which we can calculate well enough...

    return(result);
}
