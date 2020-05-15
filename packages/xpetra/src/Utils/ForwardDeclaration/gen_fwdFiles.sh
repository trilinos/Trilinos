#!/bin/bash

#FIXME JHU 22-Apr-2019 Until ETI is complete, don't run this script.
exit

package="Xpetra"
classListDir=../ClassList

for classList in ${classListDir}/*.classList
  do
  i=$(basename $classList .classList)

  tmpl=$i.tmpl

  for className in `cat $classList | grep -v \#`
    do
    uppercaseClassName=$(echo $className | tr '[a-z]' '[A-Z]')
    cat $tmpl | sed "s/\$TMPL_UPPERCASECLASS/$uppercaseClassName/g" | sed "s/\$TMPL_CLASS/$className/g" > $package"_"$className"_fwd.hpp"
  done

done
