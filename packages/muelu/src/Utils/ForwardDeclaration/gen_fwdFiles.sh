#!/bin/bash

classListDir=../ClassList/

for i in Non-Templated LO-GO-NO-LMO SC-LO-GO-NO-LMO SC-LO-GO SC-LO
  do

  classList=$classListDir/$i.classList
  tmpl=$i.tmpl

  for className in `cat $classList | grep -v \#`
    do
    uppercaseClassName=$(echo $className | tr '[a-z]' '[A-Z]')
    cat $tmpl | sed "s/\$TMPL_UPPERCASECLASS/$uppercaseClassName/g" | sed "s/\$TMPL_CLASS/$className/g" > "MueLu_"$className"_fwd.hpp"
  done

done
