C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDELB (A, IDELB, NAMELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, KATRIB, *)
C=======================================================================

C   --*** RDELB *** (GEN3D) Read database element blocks
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --RDELB reads the element block information from the database.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IDELB - OUT - the element block ID for each block
C   --   NUMELB - OUT - the number of elements for each block
C   --   NUMLNK - OUT - the number of nodes per element for each block
C   --   NUMATR - OUT - the number of attributes for each block
C   --   LINK - OUT - the connectivity array (4 nodes per element)
C   --   KATRIB - OUT - the dynamic memory pointer to the attribute array
C   --      (named 'ATRIB', packed)
C   --   * - return statement if end of file or read error
C   --
C   --Common Variables:
C   --   Uses NDBIN of /DBASE/
C   --   Uses NUMEL, NELBLK of /DBNUMS/

      INCLUDE 'exodusII.inc'
      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbnums.blk'

      DIMENSION A(*)
      CHARACTER*32 NAMELB(NELBLK)
      INTEGER IDELB(NELBLK)
      INTEGER NUMELB(NELBLK)
      INTEGER NUMLNK(NELBLK)
      INTEGER NUMATR(NELBLK)
      INTEGER LINK(4,NUMEL)

      CHARACTER*5 STRA

      IEATR = 0
      ISATR = 0
      CALL MDRSRV ('ATRIB', KATRIB, IEATR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 30

C ... Set all entries in the LINK array to -1.
C     This is a kludge, but lets us discriminate between quads, tris,
C     and shells in other parts of the program (newess)
      CALL INIINT(4*NUMEL, -1, LINK)

      call exgebi (ndbin, idelb, ierr)
      MAXNQ = 0
      DO  5 IELB = 1, NELBLK
        call exgelb(ndbin, idelb(ielb), namelb(ielb), numelb(ielb),
     &    numlnk(ielb), numatr(ielb), ierr)
        if (ierr .ne. 0) goto 30

C ... See if there are any non-quad element blocks since we have to
C     read in the connectivity different for exodusII than the code
C     expects
        if (numlnk(ielb) .ne. 4) then
          if (numelb(ielb)*numlnk(ielb) .gt. MAXNQ)
     *      MAXNQ = numelb(ielb)*numlnk(ielb)
        end if
 5    continue

      if (maxnq .gt. 0) then
        call mdrsrv('LINTMP', klntmp, maxnq)
        call mdstat (nerr, mem)
        IF (NERR .GT. 0) GOTO 20
      end if

      IEND = 0
      DO 10 IELB = 1, NELBLK
         ISTART = IEND + 1
         IEND = ISTART + NUMELB(IELB) - 1
         IF (IEND .GT. NUMEL) GOTO 40

         if (numatr(ielb) .gt. 0) then
            ISATR = IEATR + 1
            IEATR = ISATR + NUMATR(IELB) * NUMELB(IELB) - 1
            CALL MDLONG ('ATRIB', KATRIB, IEATR)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 20
         end if
         CALL RELBLK (IDELB(IELB), NUMELB(IELB), NUMLNK(IELB),
     &        NUMATR(IELB), LINK(1,ISTART), A(KATRIB+ISATR-1),
     &        a(klntmp), *50)
   10 CONTINUE

   20 CONTINUE
      if (maxnq .gt. 0) then
        call mddel('LINTMP')
        call mdstat (nerr, mem)
        IF (NERR .GT. 0) GOTO 50
      end if
      RETURN

   30 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENT BLOCK SIZING PARAMETERS for block '
     &   // STRA(:LSTRA))
      GOTO 50
   40 CONTINUE
      CALL PRTERR ('FATAL', 'Number of elements in blocks > NUMEL')
      GOTO 50
   50 CONTINUE
      RETURN 1
      END
