#!/bin/sh
for i in `ls *.c *.h`; do
    echo $i
    clang-format-mp-3.9 -style=file $i > tmp
    mv tmp $i
done
