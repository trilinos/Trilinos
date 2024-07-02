
#-----------------------------------------------------------------------------
#  Hardware locality detection and control library.
#
#  Acquisition information:
#    Date checked:  November 2011
#    Checked by:    H. Carter Edwards <hcedwar AT sandia.gov>
#    Source:        http://www.open-mpi.org/projects/hwloc/
#    Version:       1.3
#

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HWLOC
  REQUIRED_HEADERS hwloc.h
  REQUIRED_LIBS_NAMES "hwloc"
  )

