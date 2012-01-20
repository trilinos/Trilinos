#!/bin/bash

# test if fwd headers compile.

headers=$(find ../../ForwardDeclaration/ -name "*.hpp" -exec basename {} \;)

IFS=$'\n'
for header in $headers; do
    baseName=$(basename $header .hpp)
    file=$(echo $baseName".cpp")

    baseName2=$(basename $header _fwd.hpp)
    header2=$(echo $baseName2".hpp")

    echo "#include \"$header\""           >   $file
    echo "#include \"$header2\""          >>  $file
done
