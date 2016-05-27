# dev

* converted to use roxygen2

* Added a `NEWS.md` file to track changes to the package.

# bsam 0.43.1

* ported from source 2016-05-27 mdsumner@gmail.com

```R
f <- "http://web.science.mq.edu.au/~ijonsen/code/bsam_0.43.1.tar.gz"
download.file(f, basename(f), mode = "wb")
system(sprintf("tar zxvf %s", basename(f)))
Rd2roxygen::Rd2roxygen("bsam")
```

