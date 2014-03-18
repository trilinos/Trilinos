#!/usr/bin/env python
#ML timelines for matrix matrix multiply in prolongator smoothing
#*note* removelines.sh must be in a directory in your path.  arcadia changes directories, so "PATH=blah:." may not be enough.
LABELS    = ['i&x', 'serialCore', 'fc']
TIMELINES    = ['AB pre-multiply', 'AB multiply time', 'AB post-multiply']
def PARSEFUNC(fileName,stringToFind):
  # ML timings for all two matrix mults except those associated with repartitioning
  return "removeLines.sh -p \"product AB .Pmat_|product AB .rebalance|product AB .Rebal\" -n 5 " + fileName + " | grep -i \"" + stringToFind + "\" | cut -f2 -d'='"
