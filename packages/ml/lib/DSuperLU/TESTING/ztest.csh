#!/bin/csh

set ofile = ztest.out			# output file
if ( -e $ofile ) then
    rm -f $ofile
endif
echo "Double-precision complex testing output" > $ofile


set MATRICES     = (LAPACK)
set NVAL         = (9 19)
set NRHS         = (5)
set LWORK        = (0 10000000)
set PANELSIZE    = (8)
set RELAX        = (4)

#
# Loop through all matrices ...
#
foreach m ($MATRICES)

  #--------------------------------------------
  # Test matrix types generated in LAPACK-style
  #--------------------------------------------
  if  ($m == 'LAPACK') then
      echo '== LAPACK test matrices' >> $ofile
      foreach n ($NVAL)
        foreach s ($NRHS)
          foreach l ($LWORK)
	    echo '' >> $ofile
            echo 'n='$n 'nrhs='$s 'lwork='$l >> $ofile
            sp_lintstz -t "LA" -l $l -n $n -s $s >> $ofile
          end
        end
      end
  #--------------------------------------------
  # Test a specified sparse matrix
  #--------------------------------------------
  else
    echo '' >> $ofile
    echo '== sparse matrix:' $m >> $ofile
    foreach  w ($PANELSIZE)
      foreach r ($RELAX)
        foreach s ($NRHS)
          foreach l ($LWORK)
	    echo '' >> $ofile
            echo 'w='$w 'relax='$r 'nrhs='$s 'lwork='$l >> $ofile
            sp_lintstz -t "SP" -w $w -r $r -s $s -l $l \
                             < ../EXAMPLE/$m >> $ofile
	  end
        end
      end
    end
  endif

end


