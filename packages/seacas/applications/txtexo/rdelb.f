C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDELB (NTXT, IELB, IDELB, NUMELB, NUMLNK, NUMATR,
     &   NAMELB, A, KLINK, KATRIB, *)
C=======================================================================

C   --*** RDELB *** (TXTEXO) Read database element block
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDELB reads the element block information from the text file.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   IELB - IN - the element block number
C   --   IDELB - IN - the id for this block
C   --   NUMELB - IN - the number of elements in this block
C   --   NUMLNK - IN - the number of nodes per element in this block
C   --   NUMATR - IN - the number of attributes in this block
C   --   A - IN - the dynamic memory base array
C   --   KLINK - IN - pointer to the element connectivity for this block
C   --   KATRIB - IN - pointer to the attributes for this block
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      include 'exodusII.inc'
      DIMENSION A(*)
      CHARACTER*(MXSTLN) NAMELB
      CHARACTER*32 STRA

      NAMELB = ' '
      READ (NTXT, *, END=110, ERR=110)
      READ (NTXT, *, END=110, ERR=110) IDELB, NUMELB, NAMELB

C ... Strip everything in namelb from first space to end
      IEX = index(namelb, " ")
      if (iex .gt. 0) then
        namelb(iex:) = " "
      end if

      READ (NTXT, *, END=110, ERR=110) NUMLNK, NUMATR

      IECON = NUMLNK * NUMELB
      CALL MDLONG ('LINK', KLINK, IECON)
      IEATR = NUMATR * NUMELB
      CALL MDLONG ('ATRIB', KATRIB, IEATR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      CALL RDEB1 (NTXT, IELB, NUMELB, NUMLNK, NUMATR,
     &   A(KLINK), A(KATRIB), max(1,numatr), *120)

  100 CONTINUE
      RETURN

  110 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENT BLOCK SIZING PARAMETERS for block '
     &   // STRA(:LSTRA))
  120 CONTINUE
      RETURN 1
      END
