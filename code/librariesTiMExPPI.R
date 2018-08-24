pkgLoad <- function(x)
{
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    if (!require(x, character.only = TRUE))
      stop("Package not found")
  }
}

pkgLoad("TiMEx")
pkgLoad("igraph")
pkgLoad("cplexAPI")
pkgLoad("gtools")
