      program testcpd

c
c This is a test program for the Fortran binding of the EXODUS II
c database copy function (excopy).
c
      implicit none

      include 'exodusII.inc'

      integer iin, iout, exoid, exoid1, ierr, cpu_ws, io_ws

      real vers

      data iin /5/, iout /6/

c
c open EXODUS II input file
c

c the setting of cpu_ws isn't used for copying but will test the
c conversion routines

      cpu_ws = 8
      io_ws = 4

      exoid = exopen ("test.exo", EXREAD, cpu_ws, io_ws, vers, ierr)
      write (iout, '(/"after exopen, error = ",i3)')
     1			ierr

      write (iout, '("test.exo is an EXODUSII file; version ",
     1                f4.2)') vers
      write (iout, '(" I/O word size: ",i4)') io_ws

c
c  create EXODUS II output file with default size reals
c
c the setting of cpu_ws isn't used for copying but will test the
c conversion routines

      cpu_ws = 8
      io_ws = 0

      exoid1 = excre ("testcp.exo",
     1               EXCLOB, cpu_ws, io_ws, ierr)
      write (iout,'("after excre, id = ", i3, ", error = ",i3)') 
     1               exoid1, ierr
      write (iout,'(" I/O word size: ",i4)') io_ws

      write (iout,'("after excre, error = ", i4)') ierr

      call excopy (exoid, exoid1, ierr)
      write (iout, '(/"after excopy, error = ", i3)' ) ierr

      call exclos (exoid, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr

      call exclos (exoid1, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr

      stop
      end

