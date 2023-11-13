#!/bin/sh
for i in $(ls ./*.[cCh]); do
    echo $i
    /usr/local/bin/clang-format -i -style=file $i
done

