gviz_rescale <- function(z, ncolor=100, xrange=range(z, na.rm=TRUE))  #https://github.com/jotsetung/Gviz/blob/master/R/Gviz-methods.R  -> function .z2icol
{
  res <- round((z - xrange[1])/diff(xrange) * (ncolor - 1)) + 1
  res[res > ncolor] = ncolor
  res[res < 1] = 1
  return(res)
}