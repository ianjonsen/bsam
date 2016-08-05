# dev

# bsam 0.43.1 (pre-CRAN release)

* ported from source 2016-05-27 mdsumner@gmail.com

* converted to use roxygen2

* Added a `NEWS.md` file to track changes to the package.

```R
f <- "http://web.science.mq.edu.au/~ijonsen/code/bsam_0.43.1.tar.gz"
download.file(f, basename(f), mode = "wb")
system(sprintf("tar zxvf %s", basename(f)))
Rd2roxygen::Rd2roxygen("bsam")
```

# 1.0.0 

* 





