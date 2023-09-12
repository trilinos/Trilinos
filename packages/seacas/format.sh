#!/bin/sh
for i in $(ls ./*.[cCh]); do
    echo $i
    /opt/homebrew/bin/clang-format -i -style=file $i
done

