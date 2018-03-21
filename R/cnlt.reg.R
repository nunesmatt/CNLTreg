cnlt.reg <-
function (x, f, P = 50, returnall = FALSE, nkeep=2, ...) 
{
    n <- length(x)
    vec <- matrix(0, P, n - nkeep)
    deni <- df <- NULL
    aveghat <- matrix(0,1,n)

    for (i in 1:P) {
        cat(i, "...\n")
        v <- sample(1:n, (n - nkeep), FALSE)
        vec[i, ] <- v
        deni <- denoisepermC(x, f, nkeep=nkeep,mod = v, returnall = FALSE, ...)
        aveghat <- aveghat + matrix(deni,nrow=1)
    }
    aveghat <- aveghat/P
    if (returnall) {
        return(list(vec = vec, aveghat = aveghat))
    }
    else {
        return(aveghat)
    }
}
