# bsam


NOTE: not yet working or supported

- work in progress sourced from http://web.science.mq.edu.au/~ijonsen/code.html
- convert to roxygen2

todo
- set up for devtools install
- pass check
- add tests
- reduce example run time
- "tod" param was in Rd but not in dat4jags
- simTrack has undocumented params, including "T"
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
  

checking R code for possible problems ... NOTE
Found an obsolete/platform-specific call in the following functions:
  'diagSSM' 'plotSSM'
Found the platform-specific devices:
  'quartz' 'windows'
dev.new() is the preferred way to open a new device, in the unlikely
event one is needed.

