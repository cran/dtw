###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>               #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: triangleFixing.R 437 2018-05-17 14:47:39Z tonig $
#                                                             #
###############################################################


########################################
## Fix a distance matrix in order to obey the triangle inequality.
## Triangle fixing algorithm - implementation of algorithm 3.1
## (Metric_Nearness_L2) in Brickell, J., Dhillon, I., Sra, S., and
## Tropp, J. (2008). The Metric Nearness Problem. SIAM. J. Matrix
## Anal. & Appl. 30, 375-396.
##
## Wrapper to the native function


`triangleFixing` <-
function(D, tolerance=1e-15, max.iter=1e6) {

    ## D may be a lower-triangular dist object
    Ds <- as.matrix(D)                  
    n <- nrow(Ds)

    out <- .C("triangle_fixing_l2",
            ## IN+OUT
            M=as.double(Ds),
            iter=as.integer(max.iter),
            ## IN
            as.integer(n),
            as.double(tolerance),
            ## OUT
            delta=as.double(0.0))

    if(out$iter==0)
        warning("Convergence not reached")

    # cat(max.iter-out$iter)
    
    dim(out$M) <- c(n,n)
    return(out$M)
}


## Adapted from the fossil package
tri.ineq.show <- function (D) 
{
    mat <- as.matrix(D)
    n <- dim(mat)[1]
    ineq.idx <- c()
    for (i in 1:(n - 2)) {
        for (j in (i + 1):(n - 1)) {
            for (k in (j + 1):n) {
                sds <- c(mat[j, i], mat[k, i], mat[k, j])
                lng <- max(sds)
                if (lng > (sum(sds) - lng)) {
                    # cat(sprintf("Check triangle %d %d %d\n",i,j,k))
                    ineq.idx <- rbind(ineq.idx,c(i,j,k))
                }
            }
        }
    }
    return(ineq.idx)
}
