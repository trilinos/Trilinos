
#-----------------------------------------------------------------------------
#  Hardware locality detection and control library.
#
#  Acquisition information:
#    Date checked:  July 2014
#    Checked by:    H. Carter Edwards <hcedwar AT sandia.gov>
#    Source:        https://code.google.com/p/qthreads
#

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( QTHREAD
  REQUIRED_HEADERS qthread.h
  REQUIRED_LIBS_NAMES "qthread"
  )

