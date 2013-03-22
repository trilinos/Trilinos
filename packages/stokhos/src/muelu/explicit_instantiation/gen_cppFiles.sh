#!/bin/bash

IFS=$'\n'

classListDir=../../../../muelu/src/Utils/ClassList/

for i in SC-LO-GO-NO-LMO SC-LO-GO SC-LO
  do

  classList=$classListDir/$i.classList
  tmpl=$i.tmpl

  for className in `cat $classList | grep -v ^\# | cut -d "-" -f1 | sed 's/ //'`
    do
    if ! grep -q $className $classListDir/EI-Exceptions.classList && ! grep -q $className exceptions.classList
        then

        condition=$(cat $classList | grep "^$className -" | cut -d "-" -f2-)

        if [ -n "$condition" ]; then
            conditionOpen1=$(echo "#include \\\"MueLu_ConfigDefs.hpp\\\"")
            conditionOpen2=$(echo $condition | sed 's/^[ ]*//' | sed 's/\&/\\\&/g')
            conditionClose="#endif"
        else
            conditionOpen1=""
            conditionOpen2=""
            conditionClose=""
        fi

        cat $tmpl \
            | sed "s/\$TMPL_CLASS/$className/g" \
            | sed "s/\$TMPL_CONDITION_OPEN1/$conditionOpen1/g" \
            | sed "s/\$TMPL_CONDITION_OPEN2/$conditionOpen2/g" \
            | sed "s/\$TMPL_CONDITION_CLOSE/$conditionClose/g" \
            > MueLu_$className.cpp
    fi

  done

done
