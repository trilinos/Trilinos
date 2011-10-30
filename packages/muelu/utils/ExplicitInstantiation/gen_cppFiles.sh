#!/bin/bash

for i in "LO-GO-NO-LMO"
  do

  classList=$i.classList
  tmpl=$i.tmpl
  
  for className in `cat $classList`
    do
    cat $tmpl | sed "s/\$TMPL_CLASS/$className/g" > $className.cpp
  done

done