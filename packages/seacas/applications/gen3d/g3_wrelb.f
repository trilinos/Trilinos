C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRELB (A, IA, BLKTYP, NAMELB, IBPARM,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB, ELATTR,
     &   IXEL, INCEL, NREL, IELCOL, IXNP, NRNP, NUMELB3)
C=======================================================================

C   --*** WRELB *** (GEN3D) Write 3D element blocks
C   --   Written by Amy Gilkey - revised 09/02/87
C   --
C   --WRELB calculates and writes the element block information for the
C   --3D database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   BLKTYP - IN - the element block type
C   --   IBPARM - IN - the block parameters (defined by the block type)
C   --   IDELB - IN - the id for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the 2D connectivity array (always 4 nodes)
C   --   ATRIB - IN - the 2D attribute array (packed)
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses NUMEL3 of /DBNUM3/
C   --   Uses NEREPL of /PARAMS/

      INCLUDE 'exodusII.inc'
      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xyzmir.blk'

      REAL A(*)
      INTEGER IA(*)
      CHARACTER BLKTYP(NELBLK)
      CHARACTER*(MXSTLN) NAMELB(NELBLK)
      INTEGER IBPARM(4,NELBLK)
      INTEGER IDELB(NELBLK)
      INTEGER NUMELB(NELBLK)
      INTEGER NUMELB3(NELBLK)
      INTEGER NUMLNK(NELBLK)
      INTEGER NUMATR(NELBLK)
      INTEGER LINK(4,*)
      REAL ATRIB(*)
      REAL ELATTR(*)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)
      INTEGER IXNP(*), NRNP(*)

C ... Update block names
      DO I = 1, NELBLK
        IF (NAMELB(I)(:4) .EQ. 'QUAD' .OR.
     *    NAMELB(I)(:4) .EQ. 'quad' .OR.
     *    NAMELB(I) .EQ. ' ') then
          NAMELB(I) = 'HEX'
        else IF (NAMELB(I)(:3) .EQ. 'BAR'  .OR.
     *      NAMELB(I)(:3) .EQ. 'bar'  .OR.
     *      NAMELB(I)(:5) .EQ. 'TRUSS'.OR.
     *      NAMELB(I)(:5) .EQ. 'truss'.OR.
     *      NAMELB(I)(:4) .EQ. 'BEAM' .OR.
     *      NAMELB(I)(:4) .EQ. 'beam') then
          NAMELB(I) = 'SHELL'
        else if (NAMELB(I)(:3) .EQ. 'TRI') then
          NAMELB(I) = 'WEDGE'
        end if
      end do

      MAXID = 0
      MAXEL = 0
      MAXATR = 0
      DO IELB = 1, NELBLK
        MAXID = MAX (MAXID, IDELB(IELB))
        MAXEL = MAX (MAXEL, NUMELB(IELB))
        MAXATR = MAX (MAXATR, NUMATR(IELB))
      end do
      CALL MDRSRV ('LINK3', KLINK3, 8 * MAXEL*NEREPL)
      CALL MDRSRV ('ATRIB3', KATR3, MAXATR * MAXEL*NEREPL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 50

      IEL = 1
      IATR = 1

      DO IELB = 1, NELBLK

C      --Block id and start/end for blocks that become multiple blocks
         NRSTRT = 1
         IDEBL3 = IDELB(IELB)
         IBLK = 1
         LSTID = 1
         IF (BLKTYP(IELB) .EQ. 'T') THEN
            NREND = IBPARM(1,IELB) - 1
            NBLK = IBPARM(4,IELB)
         ELSE IF (BLKTYP(IELB) .EQ. 'S') THEN
            NREND = NRTRAN(1)
            NBLK = IBPARM(4,IELB)
         ELSE
            NREND = NEREPL
            NBLK = 1
         END IF
         IF (NBLK .GT. 1) THEN
            LSTID = MAXID
            MAXID = MAXID + NBLK-1
         END IF

   20    CONTINUE

C      --Number of elements in the block

         IF (BLKTYP(IELB) .EQ. 'C') THEN
           NUMELB3(IELB) = 0
           DO  I = 1, NUMELB(IELB)
             NUMELB3(IELB) = NUMELB3(IELB) + NREL(IEL+I-1)
           end do
         ELSE
           NUMELB3(IELB) = NUMELB(IELB) * (NREND - NRSTRT + 1)
         END IF

C      --Number of nodes per element

         NUMLN3 = NUMLNK(IELB) * 2

C      --Connectivity and attributes

C     --Change attributes of input block if they were changed by user
         IF (ELATTR(IELB) .ne. 0.0) then
            call INIREA (NUMELB(IELB) * NUMATR(IELB), ELATTR(IELB),
     &           ATRIB(IATR))
         END IF

         CALL NEWEL1 (BLKTYP(IELB),
     &        NUMELB(IELB), NUMELB3(IELB), NRSTRT, NREND,
     &        NUMLNK(IELB), NUMLN3, NUMATR(IELB), MAX(NUMATR(IELB),1),
     $        NUMNP, NUMNP3, LINK(1,IEL), IA(KLINK3), ATRIB(IATR),
     $        A(KATR3), IXEL(IEL), INCEL(IEL),
     &        NREL(IEL), IELCOL(IEL),
     $        IXNP, NRNP)

C      --Fixup connectivity if mirrored

         IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
            CALL DBMIRR (IELB, IELB, IDEBL3, NUMELB3(IELB), NUMLN3,
     *       IA(KLINK3))
         END IF

C      --Write 3D
         call expelb(ndbout, IDEBL3, NAMELB(IELB), NUMELB3(IELB),
     *     NUMLN3, NUMATR(IELB), ierr)
         if (ierr .lt. 0) then
            call exerr('gen3d2', 'Error from expelb', exlmsg)
            go to 50
         endif
         call expelc(ndbout, IDEBL3, IA(KLINK3), ierr)
         if (ierr .lt. 0) then
            call exerr('gen3d2', 'Error from expelc', exlmsg)
            go to 50
         endif
         if (numatr(ielb) .gt. 0) then
            call expeat(ndbout, IDEBL3, A(KATR3), ierr)
            if (ierr .lt. 0) then
               call exerr('gen3d2', 'Error from expeat', exlmsg)
               go to 50
            endif
         end if

         IF (IBLK .LT. NBLK) THEN
            IBLK = IBLK + 1
            LSTID = LSTID + 1
            IDEBL3 = LSTID
            NRSTRT = NREND + 1
            IF (BLKTYP(IELB) .EQ. 'T') THEN
               IF (NRSTRT .LE. IBPARM(2,IELB)) THEN
                  NREND = MIN (IBPARM(2,IELB), NREND + IBPARM(3,IELB))
               ELSE
                  NREND = NEREPL
               END IF
            ELSE IF (BLKTYP(IELB) .EQ. 'S') THEN
               NREND = NRSTRT + NRTRAN(IBLK) - 1
            END IF
            GOTO 20
         END IF

         IEL = IEL + NUMELB(IELB)
         IATR = IATR + NUMATR(IELB) * NUMELB(IELB)
       end do

   50 CONTINUE
      CALL MDDEL ('LINK3')
      CALL MDDEL ('ATRIB3')

      RETURN
      END
