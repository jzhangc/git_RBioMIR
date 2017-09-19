.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
                        Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.

                        (test only) For novel miRNA discovery, be sure to copy training.smir and test.smir files
                        to your working directory, and install pcregrep and parallel commands.")
  return(TRUE)
}
