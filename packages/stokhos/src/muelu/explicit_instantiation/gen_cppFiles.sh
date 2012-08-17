#!/bin/bash

for i in SC-LO-GO-NO-LMO SC-LO-GO SC-LO
do

  classList=../../../../muelu/src/Utils/ExplicitInstantiation/$i.classList
  tmpl=$i.tmpl
  
  for className in `cat $classList | grep -v \#`
  do
    if ! grep -q $className exceptions.classList
    then
        cat $tmpl | sed "s/\$TMPL_CLASS/$className/g" > MueLu_$className.cpp
    fi
  done

done
