#!/bin/bash

#
# Ordinal
#

classListDir=../Utils/ClassList/

echo "// Type definitions for templated classes (generally graph-related) that do not require a scalar." > MueLu_UseShortNamesOrdinal.hpp
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

echo "// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node of the current context." > MueLu_UseShortNamesScalar.hpp
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

# Add the matlab utilities to end of file
echo "#ifdef MUELU_TWOLEVELMATLABFACTORY_SHORT" >> MueLu_UseShortNamesScalar.hpp
echo "typedef MueLu::TwoLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TwoLevelMatlabFactory;" >> MueLu_UseShortNamesScalar.hpp
echo "#endif" >> MueLu_UseShortNamesScalar.hpp
echo "#ifdef MUELU_SINGLELEVELMATLABFACTORY_SHORT" >> MueLu_UseShortNamesScalar.hpp
echo "typedef MueLu::SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> SingleLevelMatlabFactory;" >> MueLu_UseShortNamesScalar.hpp
echo "#endif" >> MueLu_UseShortNamesScalar.hpp
echo "#ifdef MUELU_MATLABSMOOTHER_SHORT" >> MueLu_UseShortNamesScalar.hpp
echo "typedef MueLu::MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatlabSmoother;" >> MueLu_UseShortNamesScalar.hpp
echo "#endif" >> MueLu_UseShortNamesScalar.hpp
