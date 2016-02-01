#!/bin/bash

IFS=$'\n'

classListDir=../ClassList/

rm -Rf ETI_LO_GO_NO_classes.cmake
rm -Rf ETI_SC_LO_GO_NO_classes.cmake

  classList=$classListDir/LO-GO-NO.classList
  #for className in `cat $classList | grep -v ^\# | cut -d "-" -f1 | sed 's/ //'`
  for className in `cat $classList | grep -v ^\# | sed 's/ //'`
    do

    if ! grep -q -x $className $classListDir/EI-Exceptions.classList
        then

        printClassName=$(echo $className | sed 's/ /./g' | sed 's/#/?/g' | sed 's/(/[/g' | sed 's/)/]/g')
	echo "APPEND_SET(MUELU_LO_GO_NO_ETI_CLASSES MueLu::$printClassName )" >> ETI_LO_GO_NO_classes.cmake
    fi

  done


  classList=$classListDir/SC-LO-GO-NO.classList
  for className in `cat $classList | grep -v ^\# | sed 's/ //'`
    do

    if ! grep -q -x $className $classListDir/EI-Exceptions.classList
        then

        printClassName=$(echo $className | sed 's/ /./g' | sed 's/#/?/g' | sed 's/(/[/g' | sed 's/)/]/g')
	echo "APPEND_SET(MUELU_SC_LO_GO_NO_ETI_CLASSES MueLu::$printClassName )" \
		>> ETI_SC_LO_GO_NO_classes.cmake
    fi

  done

