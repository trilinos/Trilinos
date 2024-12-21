#!/bin/bash

#
# Ordinal
#

classListDir=../Utils/ClassList/

echo "// @HEADER" > MueLu_UseShortNamesOrdinal.hpp
echo "// *****************************************************************************" >> MueLu_UseShortNamesOrdinal.hpp
echo "//        MueLu: A package for multigrid based preconditioning" >> MueLu_UseShortNamesOrdinal.hpp
echo "//" >> MueLu_UseShortNamesOrdinal.hpp
echo "// Copyright 2012 NTESS and the MueLu contributors." >> MueLu_UseShortNamesOrdinal.hpp
echo "// SPDX-License-Identifier: BSD-3-Clause" >> MueLu_UseShortNamesOrdinal.hpp
echo "// *****************************************************************************" >> MueLu_UseShortNamesOrdinal.hpp
echo "// @HEADER" >> MueLu_UseShortNamesOrdinal.hpp
echo "" >> MueLu_UseShortNamesOrdinal.hpp
echo "// Type definitions for templated classes (generally graph-related) that do not require a scalar." >> MueLu_UseShortNamesOrdinal.hpp
echo >> MueLu_UseShortNamesOrdinal.hpp
echo "#include <Xpetra_UseShortNamesOrdinal.hpp>" >> MueLu_UseShortNamesOrdinal.hpp
echo >> MueLu_UseShortNamesOrdinal.hpp

for i in LO-GO-NO Non-Templated
  do

  classList=$classListDir/$i.classList
  tmpl=$i.tmpl

  for className in `cat $classList | grep -v ^\# | cut -d "-" -f1 | sed 's/ //'`
    do
    uppercaseClassName=$(echo $className | tr '[a-z]' '[A-Z]')
    cat $tmpl | sed "s/\$TMPL_UPPERCASECLASS/$uppercaseClassName/g" | sed "s/\$TMPL_CLASS/$className/g" >> MueLu_UseShortNamesOrdinal.hpp
  done
done

#
# Scalar
#
echo "// @HEADER" > MueLu_UseShortNamesScalar.hpp
echo "// *****************************************************************************" >> MueLu_UseShortNamesScalar.hpp
echo "//        MueLu: A package for multigrid based preconditioning" >> MueLu_UseShortNamesScalar.hpp
echo "//" >> MueLu_UseShortNamesScalar.hpp
echo "// Copyright 2012 NTESS and the MueLu contributors." >> MueLu_UseShortNamesScalar.hpp
echo "// SPDX-License-Identifier: BSD-3-Clause" >> MueLu_UseShortNamesScalar.hpp
echo "// *****************************************************************************" >> MueLu_UseShortNamesScalar.hpp
echo "// @HEADER" >> MueLu_UseShortNamesScalar.hpp
echo "" >> MueLu_UseShortNamesScalar.hpp
echo "// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node of the current context." >> MueLu_UseShortNamesScalar.hpp
echo >> MueLu_UseShortNamesScalar.hpp
echo "#include <Xpetra_UseShortNamesScalar.hpp>" >> MueLu_UseShortNamesScalar.hpp
echo >> MueLu_UseShortNamesScalar.hpp

i=SC-LO-GO-NO
classList=$classListDir/$i.classList
tmpl=$i.tmpl

for className in `cat $classList | grep -v ^\# | cut -d "-" -f1 | sed 's/ //'`
  do
  uppercaseClassName=$(echo $className | tr '[a-z]' '[A-Z]')
  cat $tmpl | sed "s/\$TMPL_UPPERCASECLASS/$uppercaseClassName/g" | sed "s/\$TMPL_CLASS/$className/g" >> MueLu_UseShortNamesScalar.hpp
done

# AmesosSmoother and IfpackSmoother are special (they need only one template parameter)
echo "#ifdef MUELU_AMESOSSMOOTHER_SHORT" >> MueLu_UseShortNamesOrdinal.hpp
echo "typedef MueLu::AmesosSmoother<Node> AmesosSmoother;" >> MueLu_UseShortNamesOrdinal.hpp
echo "#endif" >> MueLu_UseShortNamesOrdinal.hpp
echo "#ifdef MUELU_IFPACKSMOOTHER_SHORT" >> MueLu_UseShortNamesOrdinal.hpp
echo "typedef MueLu::IfpackSmoother<Node> IfpackSmoother;" >> MueLu_UseShortNamesOrdinal.hpp
echo "#endif" >> MueLu_UseShortNamesOrdinal.hpp
