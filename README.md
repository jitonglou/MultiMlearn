MultiMlearn
================

The goal of `MultiMlearn` is to provide a computing approach to estimate
individualized treatment rules under the setting of multicategory
treatments. This package is designed for data from observational studies
(such as electronic health records), but it can also be used for
randomized controlled trials.

## Installation

This package is under development, so please install the latest version
from Github:

``` r
if (!requireNamespace("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jitonglou/MultiMlearn")
```

The first time installation may need to restart the R session several
times to install/update some required packages and can take up to 30
minutes. Alternatively, if you do have difficulties in the installation,
you can clone this repository to your device, and manually read the
functions in R:

``` r
files = list.files("./R", full.names = TRUE)
for (file in files){source(file)}
```

## Vignette

An implementation of the `MultiMlearn` package to a simulated dataset
can refer to [this
website](https://jitonglou.github.io/MultiMlearn/articles/example.html)
and [this R Markdown
file](https://github.com/jitonglou/MultiMlearn/blob/master/vignettes/example.Rmd).
Details of the simulation study can be found in Section S.1 of [this pdf
document](https://github.com/jitonglou/MultiMlearn/blob/master/doc/supp_v3.pdf).

## Help files of functions

The help files of functions can be found on [this
website](https://jitonglou.github.io/MultiMlearn/reference/index.html).
Also, after you successfully install the package, you can access the
help files of some main functions by:

``` r
library(MultiMlearn)
?simulate_data
?rfcv2
?mlearn.wsvm.cv
```

Otherwise, you can read the annotations in the source files:
[simdata.R](https://github.com/jitonglou/MultiMlearn/blob/master/R/simdata.R),
[rfcv2.R](https://github.com/jitonglou/MultiMlearn/blob/master/R/rfcv2.R),
and
[mlearn.R](https://github.com/jitonglou/MultiMlearn/blob/master/R/mlearn.R).
