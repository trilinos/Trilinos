C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

      program testcpnl

c
c This is a test program for the Fortran binding of the EXODUS II
c database copy function (excopy).
c
      implicit none

      include 'exodusII.inc'

      integer iin, iout, exoid, exoid1, ierr, cpu_ws, io_ws, mod_sz

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
      mod_sz = exlgmd(exoid)
      write (iout, '(" Model Size",i2)') mod_sz

c
c  create EXODUS II output file with default size reals
c
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

