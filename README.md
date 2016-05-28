# bsam


NOTE: not yet working or supported

- work in progress sourced from http://web.science.mq.edu.au/~ijonsen/code.html

```R
devtools::install_github("ijonsen/bsam")
```

# R CMD check

all passing except for quartz/windows usage - probabably can remove these tests and just use dev.new

checking R code for possible problems ... NOTE
Found an obsolete/platform-specific call in the following functions:
  'diagSSM' 'plotSSM'
Found the platform-specific devices:
  'quartz' 'windows'
dev.new() is the preferred way to open a new device, in the unlikely
event one is needed.


# todo
- add tests
- reduce example run time
- use Authors@R
- document lbt properly

## function argument issues

- "tod" param was in Rd but not in dat4jags
- simTrack has undocumented params, including "T" which maybe should be renamed
- Undocumented arguments in documentation object 'ssm'
  'adapt' 'samples' 'thin' 'chains'
  Argument items with no description in Rd object 'diagSSM':
  'fit.in' 'save.to.pdf'

Argument items with no description in Rd object 'hssm':
  'loc.list' 'model' 'adapt' 'samples' 'thin' 'chains' '...'

Argument items with no description in Rd object 'plotSSM':
  'fit.in' 'save.to.pdf'

Argument items with no description in Rd object 'simTrack':
  'theta' 'gamma' 'alpha' 'vcov'
  



