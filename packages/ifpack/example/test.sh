#!/bin/sh

echo "compiling the examples..."
make
echo "running the examples... (this may take a while)"
echo ""
FILELIST="Ifpack_ex_Amesos.exe \
          Ifpack_ex_BlockRelaxation.exe \
          Ifpack_ex_Reordering.exe \
          Ifpack_ex_Factory.exe \
          Ifpack_ex_SingletonFilter.exe
          Ifpack_ex_Filtering.exe \
          Ifpack_ex_ICT.exe"

status="true"
for i in $FILELIST
do
  echo -n "testing $i..."
  ./$i 2>&1 > /dev/null
  if [[ $? -eq 0 ]]; then
     echo " passed!"
  else
    echo " FAILED! Please check the documentation of"
    echo "IFPACK, and the configuration parameters"
    echo ""
    status="false"
  fi
done

if [[ $status = "true" ]]; then
  echo ""
  echo "All IFPACK tests passed.";
  echo ""
fi
