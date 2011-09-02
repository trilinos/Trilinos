C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
*DECK,OPNFIL
      SUBROUTINE OPNFIL
C
C     ******************************************************************
C
C     SUBROUTINE TO OPEN REQUIRED FILES
C
C     Calls subroutine ERROR
C
C     Called by MAPVAR
C
C     ******************************************************************
C
      CHARACTER*256 fntpo, fntp2, fntp3, fntp4, filnam, option, errmsg
      character*1   cdum
C
      include 'ex2tp.blk'
      include 'ntpdat.blk'
      include 'tapes.blk'
      include 'exodusII.inc'
      include 'f2kcli.inc'
C
C     ******************************************************************
C
C .. Get filename from command line.  If not specified, emit error message
      NARG = COMMAND_ARGUMENT_COUNT()
      if (narg .eq. 0) then
        CALL PRTERR ('FATAL', 'Invalid syntax.')
        CALL PRTERR ('CMDSPEC',
     *    'Syntax is: "mapvar [-output ofile.o] [-plot result_file.e]'//
     *    ' [-mesh mesh_file.g]')
        CALL PRTERR ('CMDSPEC',
     *    '                   [-interpolated int_file.int] base"')
        stop 'Syntax Error'
      end if

C ... Get 'base' and do default file assignments.
      CALL GET_COMMAND_ARGUMENT(NARG,FILNAM, LFIL, ISTATUS)
      fntpo = filnam(:lfil) // '.o'
      fntp2 = filnam(:lfil) // '.e'
      fntp3 = filnam(:lfil) // '.g'
      fntp4 = filnam(:lfil) // '.int'

C ... parse arguments to see if any files set explicitly
C     All (except for base) should be of the form "-option file"
      do iarg = 1, narg-1, 2
        CALL GET_COMMAND_ARGUMENT(iarg,  OPTION, LOPT, ISTATUS)
        CALL GET_COMMAND_ARGUMENT(iarg+1,FILNAM, LFIL, ISTATUS)
        
        if (option(1:1) .ne. '-') then
          errmsg = 'Option "'//OPTION(:LOPT)//
     *              '" does not start with "-"'
          call prterr('FATAL', errmsg(:lenstr(errmsg)))
          CALL PRTERR ('CMDSPEC',
     *'Syntax is: "mapvar [-output ofile.o] [-plot result_file.e]'//
     *' [-mesh mesh_file.g]')
          CALL PRTERR ('CMDSPEC',
     *'                   [-interpolated int_file.int] base"')
          stop 'Syntax Error'
        end if

        if (option .eq. "-output") then
          fntpo = filnam(:lfil)
        else if (option .eq. "-plot") then
          fntp2 = filnam(:lfil)
        else if (option .eq. "-mesh") then
          fntp3 = filnam(:lfil)
        else if (option .eq. "-interpolated") then
          fntp4 = filnam(:lfil)
        else
          errmsg = 'Unrecognized option "'//OPTION(:LOPT)//'"'
          call prterr('FATAL', errmsg(:lenstr(errmsg)))
          CALL PRTERR ('CMDSPEC',
     *'Syntax is: "mapvar [-output ofile.o] [-plot result_file.e]'//
     *' [-mesh mesh_file.g]')
          CALL PRTERR ('CMDSPEC',
     *'                   [-interpolated int_file.int] base"')
          stop 'Syntax Error'
        end if
      end do
      
C     OPENING OF INPUT/OUTPUT, SCRATCH AND DATA FILES
C
C     TEXT OUTPUT FILE
C
      IFILES(1)=1
      IUNIT=NTPOUT
      OPEN (UNIT=NTPOUT, FILE=fntpo(:lenstr(fntpo)), STATUS='unknown', 
     &      FORM='formatted', ERR=10)
C
C     EXODUS DATA FILE - MESH-A (MESH & SOLUTION)
C
      IFILES(3)=1
      IUNIT=NTP2
      icpuws = 0
      iows1 = 0
      ntp2ex = exopen(fntp2(:lenstr(fntp2)),EXREAD,icpuws,iows1,
     *                vers,ierr)
      if (ierr .lt. 0)then
        WRITE(NOUT,11)IERR,NTP2,fntp2(:lenstr(fntp2)),NTP2EX
        WRITE(NTPOUT,11)IERR,NTP2,fntp2(:lenstr(fntp2)),NTP2EX
        GO TO 10
      end if
 11   format(' In opnfil - error opening mesh + solution file',/,
     &' error number ',i5,' tape number ',i5,' filename ',a20,/,
     &' ntp2ex = ',i5)
C
C     GENESIS DATA FILE - MESH-B (MESH)
C
      IFILES(4)=1
      IUNIT=NTP3
      icpuws = 0
      iows2 = 0
      ntp3ex = exopen(fntp3(:lenstr(fntp3)),EXREAD,icpuws,iows2,
     *  vers,ierr)
      if (ierr .lt. 0)then
        WRITE(NOUT,12)IERR,NTP3,fntp3(:lenstr(fntp3)),NTP3EX
        WRITE(NTPOUT,12)IERR,NTP3,fntp3(:lenstr(fntp3)),NTP3EX
        GO TO 10
      end if
 12   format(' In opnfil - error opening genesis file',/,
     &' error number ',i5,' tape number ',i5,' filename ',a20,/,
     &' ntp3ex = ',i5)

C ... Find the largest name size...
      call exinq(ntp2ex, EXDBMXUSNM, namlen2, rdum, cdum, ierr)
      call exinq(ntp3ex, EXDBMXUSNM, namlen3, rdum, cdum, ierr)

      namlen = max(namlen2, namlen3)

C
C     EXODUS DATA FILE - MESH-C (MESH & INTERPOLATED SOLUTION)
C
      IFILES(5)=1
      IUNIT=NTP4
C ... Set iows for created file to default floating point word size
C     on this machine (as returned by previous call to exopen)      
      iows3 = icpuws
      icpuws = 0
      ntp4ex = excre(fntp4(:lenstr(fntp4)),EXCLOB,icpuws,
     *  iows3,ierr)
      if (ierr .lt. 0)then
        WRITE(NOUT,13)IERR,NTP4,fntp4(:lenstr(fntp4)),NTP4EX
        WRITE(NTPOUT,13)IERR,NTP4,fntp4(:lenstr(fntp4)),NTP4EX
        GO TO 10
      end if
 13   format(' In opnfil - error opening interpolation file',/,
     &' error number ',i5,' tape number ',i5,' filename ',a20,/,
     &' ntp4ex = ',i5)
C
      call exmxnm(ntp2ex, namlen, ierr)
      call exmxnm(ntp3ex, namlen, ierr)
      call exmxnm(ntp4ex, namlen, ierr)

      RETURN
C
   10 CONTINUE
      CALL ERROR ('OPNFIL','ERROR OPENING FILE','UNIT NUMBER',IUNIT,    
     1' ',0,' ',' ',1)
      END
