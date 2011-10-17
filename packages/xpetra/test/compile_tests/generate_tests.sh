#!/bin/bash
headers=$(find ../.. -name "*.hpp" -exec basename {} \;)

exclude="
Xpetra_UseShortNames.hpp
Xpetra_UseShortNamesScalar.hpp
Xpetra_UseShortNamesOrdinal.hpp

Xpetra_BlockMap.hpp
Xpetra_TpetraBlockMap.hpp
Xpetra_EpetraBlockMap.hpp

Xpetra_VbrMatrix.hpp
Xpetra_TpetraVbrMatrix.hpp

Xpetra_DoxygenDocumentation.hpp
Xpetra_UnitTestHelpers.hpp

Xpetra_RowGraph.hpp
"

not_a_templated_class=$(find ../.. -name "*Epetra*.hpp" -exec basename {} \;)
not_a_templated_class="
$not_a_templated_class

Kokkos_ConfigDefs.hpp
Kokkos_DefaultKernels.hpp
Kokkos_DefaultNode.hpp
Kokkos_SerialNode.hpp
Xpetra_ConfigDefs.hpp
Xpetra_Exceptions.hpp
Xpetra_UseDefaultTypes.hpp
Xpetra_DefaultPlatform.hpp
Xpetra_Parameters.hpp
Xpetra_TpetraConfigDefs.hpp
Xpetra_Utils.hpp
"

# headers = headers - exclude
IFS=$' '
headers=$( comm -23 <(sort <(echo $headers)) <(sort <(echo $exclude)) )
IFS=$'\n'

for header in $headers; do
    baseName=$(basename $header .hpp)
    file=$(echo $baseName".cpp")

    echo "#include \"$header\""            >  $file
done


# headers = headers - not_a_templated_class
IFS=$' '
headers=$( comm -23 <(sort <(echo $headers)) <(sort <(echo $not_a_templated_class)) )
IFS=$'\n'

for header in $headers; do
    baseName=$(basename $header .hpp)
    file=$(echo $baseName".cpp")
    className=$(sed s/_/::/ <(echo $baseName))

    echo                                   >> $file
    echo "template class $className<int>;" >> $file
done
