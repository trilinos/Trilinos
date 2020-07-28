C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,OPNFIL
      SUBROUTINE MVOPNFIL

C     ******************************************************************

C     SUBROUTINE TO OPEN REQUIRED FILES

C     Calls subroutine ERROR

C     Called by MAPVAR

C     ******************************************************************

      CHARACTER*2048 filnam, option, errmsg

      include 'exodusII.inc'
      include 'ex2tp.blk'
      include 'ntpdat.blk'
      include 'tapes.blk'
C      include 'argparse.inc'
      integer  argument_count
      external argument_count

C     ******************************************************************

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
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
      CALL get_argument(NARG,FILNAM, LFIL)
      fntpo = filnam(:lfil) // '.o'
      fntp2 = filnam(:lfil) // '.e'
      fntp3 = filnam(:lfil) // '.g'
      fntp4 = filnam(:lfil) // '.int'

C ... parse arguments to see if any files set explicitly
C     All (except for base) should be of the form "-option file"
      do iarg = 1, narg-1, 2
        CALL get_argument(iarg,  OPTION, LOPT)
        CALL get_argument(iarg+1,FILNAM, LFIL)

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

C     TEXT OUTPUT FILE

      IFILES(1)=1
      IUNIT=NTPOUT
      OPEN (UNIT=NTPOUT, FILE=fntpo(:lenstr(fntpo)), STATUS='unknown',
     &      FORM='formatted', ERR=10)

C     EXODUS DATA FILE - MESH-A (MESH & SOLUTION)

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

C     GENESIS DATA FILE - MESH-B (MESH)

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
      namlen2 = exinqi(ntp2ex, EXDBMXUSNM)
      namlen3 = exinqi(ntp3ex, EXDBMXUSNM)

      namlen = max(namlen2, namlen3)

C     EXODUS DATA FILE - MESH-C (MESH & INTERPOLATED SOLUTION)

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

      write (*,*) 'MESH-A (MESH & SOLUTION):     ',
     *  fntp2(:lenstr(fntp2))
      write (*,*) 'MESH-B (MESH):                ',
     *  fntp3(:lenstr(fntp3))
      write (*,*) 'MESH-C (MESH & INTERPOLATED): ',
     *  fntp4(:lenstr(fntp4))
      write (*,*) 'Text Output:                  ',
     *  fntpo(:lenstr(fntpo))
      call exmxnm(ntp2ex, namlen, ierr)
      call exmxnm(ntp3ex, namlen, ierr)
      call exmxnm(ntp4ex, namlen, ierr)

      RETURN

   10 CONTINUE
      CALL ERROR ('OPNFIL','ERROR OPENING FILE','UNIT NUMBER',IUNIT,
     1' ',0,' ',' ',1)
      END
