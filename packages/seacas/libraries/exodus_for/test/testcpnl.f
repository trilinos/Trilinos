C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      program testcpnl

c This is a test program for the Fortran binding of the EXODUS II
c database copy function (excopy).

      implicit none

      include 'exodusII.inc'

      integer iin, iout, exoid, exoid1, ierr, cpu_ws, io_ws, mod_sz

      real vers

      data iin /5/, iout /6/

c open EXODUS II input file

c the setting of cpu_ws isn't used for copying but will test the
c conversion routines

      cpu_ws = 8
      io_ws = 4

      exoid = exopen ("test.exo", EXREAD, cpu_ws, io_ws, vers, ierr)
      write (iout, '(/"after exopen, error = ",i3)')
     1                  ierr

      write (iout, '("test.exo is an EXODUSII file; version ",
     1                f4.2)') vers
      write (iout, '(" I/O word size: ",i4)') io_ws
      mod_sz = exlgmd(exoid)
      write (iout, '(" Model Size",i2)') mod_sz

c  create EXODUS II output file with default size reals

c the setting of cpu_ws isn't used for copying but will test the
c conversion routines

      cpu_ws = 8
      io_ws = 0

      exoid1 = excre ("testcpnl.exo",
     1               EXCLOB+EXLARG, cpu_ws, io_ws, ierr)
      write (iout,'("after excre, id = ", i3, ", error = ",i3)')
     1               exoid1, ierr
      write (iout,'(" I/O word size: ",i4)') io_ws

      mod_sz = exlgmd(exoid1)
      write (iout, '(" Model Size",i2)') mod_sz

      write (iout,'("after excre, error = ", i4)') ierr

      call excopy (exoid, exoid1, ierr)
      write (iout, '(/"after excopy, error = ", i3)' ) ierr

      call exclos (exoid, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr

      call exclos (exoid1, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr

      stop
      end

