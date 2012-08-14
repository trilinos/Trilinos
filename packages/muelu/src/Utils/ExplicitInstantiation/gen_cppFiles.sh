#!/bin/bash

for i in LO-GO-NO-LMO SC-LO-GO-NO-LMO SC-LO-GO SC-LO
  do

  classList=$i.classList
  tmpl=$i.tmpl
  
  for className in `cat $classList | grep -v \#`
    do
    cat $tmpl | sed "s/\$TMPL_CLASS/$className/g" > MueLu_$className.cpp
  done

done
