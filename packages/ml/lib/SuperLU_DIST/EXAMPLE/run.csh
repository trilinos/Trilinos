#!/bin/csh
set MATRICES = (\
af23560.rua \
bramley1.rua \
fidapm11.rua \
ex11.rua \
goodwin.rua \
memplus.rua \
onetone1.rua \
wang4.rua \
)

unalias rm

set ofile = run.out
if ( -e $ofile ) then
  rm -f $ofile
endif

touch $ofile
set matdir = /usr/tmp/xiaoye/
foreach m ($MATRICES)
    echo "Matrix" $m
    echo "" >> $ofile
    echo "++++++++" $m "++++++++" >> $ofile
    pddrive $matdir/$m >> $ofile
end
