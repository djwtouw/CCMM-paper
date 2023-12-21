save_bin <- function(val, fname)
{
    writeBin(val, fname)
}


load_bin <- function(fname, n = 1) {
    val <- readBin(fname, what = "numeric", n = n)
    return(val)
}
