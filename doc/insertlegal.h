#!/bin/tcsh
mkdir -p old
mkdir -p new

foreach f ($*)
  cp $f old
  cat $TRILINOS_HOME/doc/TrilinosCopyRightNotice.h $TRILINOS_HOME/doc/TrilinosLicenseLanguage.h $f > new/$f
  diff $f new/$f
end
