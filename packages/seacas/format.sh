#!/bin/sh
for i in $(ls *.[cCh]); do
    echo $i
    clang-format -style=file $i > tmp
    mv tmp $i
done

