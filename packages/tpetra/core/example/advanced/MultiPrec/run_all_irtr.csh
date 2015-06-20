#!/bin/tcsh
foreach mat (bcsstk14.hb bcsstk17.rsa)
  echo Matrix $mat
  foreach mach (threaded.xml serial.xml)
    echo Machine $mach
    ./Tpetra_IRTR_double.exe    --matrix-file=$mat --machine-file=$mach --param-file=irtr_double.xml
    ./Tpetra_IRTR_qd.exe        --matrix-file=$mat --machine-file=$mach --param-file=irtr_qd.xml
    ./Tpetra_IRTR_qd_double.exe --matrix-file=$mat --machine-file=$mach --param-file=irtr_qd.xml
  end
end
