#!/bin/sh

echo "compiling the examples..."
make
echo ""
echo "running the examples... (this may take a while)"
echo "NOTE: this simple test supposes that you are using MPI, and that"
echo "the command to execute MPI programs in mpirun -np <num procs> <exe>"
echo "You may need to change the source of this script."
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
