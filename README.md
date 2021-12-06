# Catchment
An R package to calculate spatial access and availability metrics.

Catchments (catchment areas) are regions defined by the effective range of a catching entity (or the result of its
catching). For example, a basin might catch water, a hospital patients, a school students, or a store or
restaurant customers. Each of these entities might have any number of conceptual catchments, such as those
externally defined (like a school district), those realized (as by the origin of actual customers), and those
projected (such as a travel-time radius).

[catchment_ratio](https://uva-bi-sdad.github.io/catchment/reference/catchment_ratio.html) is the main function,
used to define catchments by a travel-cost matrix with bounds and/or a decay function, and calculate supply to
demand ratios within them. It is a generalized implementation of a range of 2- and 3-step floating catchment
area models, which is generally meant to align with those of the [access](https://access.readthedocs.io)
Python package.

## Installation
Download R from [r-project.org](https://www.r-project.org), then install the package from an R console:

```R
# install.packages('remotes')
remotes::install_github('uva-bi-sdad/catchment')
```

And load the package:
```R
library(catchment)
```
