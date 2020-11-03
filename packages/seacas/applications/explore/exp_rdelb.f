C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDELB (NDB, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   A, C, KLINK, KATRIB, KATRNM, ISEOF, EBTYPE, EBNAME, NAMLEN)
C=======================================================================

C   --*** RDELB *** (EXPLORE) Read database element blocks
C   --
C   --RDELB reads the element block information from the database.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBLK - IN - the number of element blocks to read
C   --   IDELB - OUT - the element block ID for each block
C   --   NUMELB - OUT - the number of elements for each block
C   --   NUMLNK - OUT - the number of nodes per element for each block
C   --   NUMATR - OUT - the number of attributes for each block
C   --   A - IN - the dynamic memory base array
C   --   KLINK - OUT - the dynamic memory pointer to the connectivity array
C   --      (named 'LINK')
C   --   KATRIB - OUT - the dynamic memory pointer to the attribute array
C   --      (named 'ATRIB')
C   --   ISEOF - IN/OUT - set true if end of file read

      include 'exodusII.inc'

      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(NAMLEN) EBNAME(*)

      DIMENSION A(*)
      CHARACTER*1 C(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIINT (NELBLK, 0, IDELB)
      CALL INIINT (NELBLK, 0, NUMELB)
      CALL INIINT (NELBLK, 0, NUMLNK)
      CALL INIINT (NELBLK, 0, NUMATR)
      CALL INISTR (NELBLK, ' ', EBTYPE)
      CALL INISTR (NELBLK, ' ', EBNAME)

C ... Get element block ids
      if (nelblk .gt. 0) then
        call exgebi(ndb, idelb, ierr)
        if (ierr .ne. 0) go to 120
      end if

C ... Read element block sizing parameters
      IELNK = 0
      IEATR = 0
      INATR = 0
      DO 100 IELB = 1, NELBLK
        call exgelb(ndb, idelb(ielb), ebtype(ielb), numelb(ielb),
     &       numlnk(ielb), numatr(ielb), ierr)
        if (ierr .ne. 0) go to 120

        if (ebtype(ielb) .eq. 'nsided' .or.
     *      ebtype(ielb) .eq. 'NSIDED') THEN
          IELNK = IELNK + NUMLNK(IELB)
        else
          IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)
      end if
      IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
      INATR = INATR + NUMATR(IELB)
 100  CONTINUE

      CALL MDRSRV ('LINK',  KLINK,  IELNK)
      CALL MDRSRV ('ATRIB', KATRIB, IEATR)
      call mcrsrv ('ATRNM', KATRNM, INATR*NAMLEN)
      CALL MDSTAT (NERR, MEM)
      if (nerr .gt. 0) go to 140

C ... Read element block connectivity and attributes
      ielnk = 0
      ieatr = 0
      inatr = 0
      do 110 ielb = 1, nelblk
        islnk = ielnk + 1
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
          ielnk = islnk + numlnk(ielb) - 1
        else
          ielnk = islnk + numlnk(ielb) * numelb(ielb) - 1
        end if
        isatr = ieatr + 1
        ieatr = isatr + numatr(ielb) * numelb(ielb) - 1

        CALL RDEB1 (NDB,
     &    IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB),
     &    A(KLINK+ISLNK-1), A(KATRIB+ISATR-1),
     &    C(KATRNM+NAMLEN*INATR), NAMLEN)

        inatr = inatr + numatr(ielb)
 110  CONTINUE

C ... Read element block names (if they exist)
      CALL EXGNAMS(NDB, EXEBLK, nelblk, ebname, ierr)
      RETURN

  120 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'ELEMENT BLOCK SIZING PARAMETERS for block', IELB
      GOTO 130
  130 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
  140 CONTINUE
      RETURN

10000  FORMAT (5 (A, I12))
      END
