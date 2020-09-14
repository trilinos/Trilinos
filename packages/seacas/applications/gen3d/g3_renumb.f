C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENUMB (A, BLKTYP, NUMELB, LINK, XN, YN,
     &   IXEL, INCEL, NREL, IELCOL, IXNP, NRNP, NPCEN, IELROW, IROT)
C=======================================================================

C   --*** RENUMB *** (GEN3D) Find the new element and node numbers
C   --   Written by Amy Gilkey - revised 04/22/88
C   --
C   --RENUMB calculates the new node and element numbers and the number
C   --of nodes/element generated for each node/element.  Unless the
C   --node or element is part of a center block, the number generated
C   --is simply the number of repetitions.  For a center block, the
C   --elements and nodes are ordered by row and column, and the number
C   --of repetitions is assigned from the column number.  The new node
C   --number is assigned so that the node and all the nodes it generates
C   --are grouped together.  For example, if the new node number is 20,
C   --the node it generates in the first translation would be 21.  The
C   --elements are handled the same way.
C   --
C   --The output number of nodes and elements are set.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base
C   --   BLKTYP - IN - the element block type
C   --   NUMELB - IN - the number of elements for each block
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   XN, YN - IN - the nodal coordinates
C   --   IXEL - OUT - the new index for each element
C   --   INCEL - OUT - the increment for each element, needed for blocks that
C   --      become multiple blocks
C   --   NREL - OUT - the number of new elements generated for each element
C   --   IELCOL - OUT - the row number for each element, 0 if not needed,
C   --      needed for single point center, useful for other centers
C   --   IXNP - OUT - the new index for each node
C   --   NRNP - OUT - the number of new nodes generated for each node
C   --   NPCEN - OUT - the node numbers of the center nodes by column and row
C   --   IELROW - SCRATCH - size = NUMEL
C   --   IROT - OUT - tracks connectivity rotation to update sset faces
C   --
C   --Common Variables:
C   --   Uses NUMNP, NUMEL of /DBNUMS/
C   --   Sets NUMNP3, NUMEL3 of /DBNUM3/
C   --   Uses NNREPL, NEREPL, DIM3, NUMCOL, NUMROW of /PARAMS/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      DIMENSION A(*)
      CHARACTER BLKTYP(*)
      INTEGER NUMELB(*)
      INTEGER LINK(4,*)
      REAL XN(NUMNP), YN(NUMNP)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)
      INTEGER IXNP(*), NRNP(*)
C.............NPCEN(NUMCOL,NUMROW)
      INTEGER NPCEN(*)
      INTEGER IELROW(NUMEL)
      INTEGER IROT(NUMEL)

      LOGICAL ALLCEN

C   --Most elements generate an element for each plate/slice

      DO 10 I = 1, NUMEL
         NREL(I) = NEREPL
         IELCOL(I) = 0
         IROT(I) = 1
   10 CONTINUE

C   --Most nodes generate a node for each plate/slice (+1)

      DO 20 INP = 1, NUMNP
         NRNP(INP) = NNREPL
   20 CONTINUE

      IF (NUMCOL .GT. 0) THEN
         NEROW = 0
         NUMROW = 0

C      --Save the element numbers of center block element in IELROW

         NCEN = 0
         IEL = 1
         DO 40 IELB = 1, NELBLK
            IF (BLKTYP(IELB) .EQ. 'C') THEN
               DO 30 I = 1, NUMELB(IELB)
                  NCEN = NCEN + 1
                  IELROW(NCEN) = IEL
                  IEL = IEL + 1
   30          CONTINUE
            ELSE
               IEL = IEL + NUMELB(IELB)
            END IF
   40    CONTINUE

C      --Set up the connectivity indices
         CALL IXLINK (NCEN, IELROW, LINK, XN, YN, IROT)

C      --Order elements by rows and columns

         CALL MDRSRV ('IXROW', KIXROW, NUMEL+1)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

         CALL MAKROW (LINK, XN, NCEN, NUMCOL, NEROW, IELROW, A(KIXROW))

         CALL MDLONG ('IXROW', KIXROW, NEROW+1)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

C      --Put elements into row by column array (column 1 of each row is
C      --not necessarily the center)

         CALL MDRSRV ('IELCEN', KELCEN, NUMCOL * NEROW)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

         CALL FELCEN (NUMCOL, NEROW, IELROW, A(KIXROW), A(KELCEN))

         CALL MDDEL ('IXROW')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

C      --Determine which elements are on the center column

         CALL MDRSRV ('ICOL1', KICOL1, NEROW)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

         CALL MRKCEN (LINK, XN, NUMCOL, NEROW, A(KELCEN),
     &      A(KICOL1), ALLCEN)

C      --Find the proper columns for the elements

         CALL MDRSRV ('IRBOT', KIRBOT, NEROW)
         CALL MDRSRV ('IRTOP', KIRTOP, NEROW)

         CALL MDRSRV ('IRBDIF', KIBDIF, NEROW)
         CALL MDRSRV ('NPBROW', KNPBR, NUMNP)
         CALL MDRSRV ('NPBCOL', KNPBC, NUMNP)
         CALL MDRSRV ('NPTROW', KNPTR, NUMNP)
         CALL MDRSRV ('NPTCOL', KNPTC, NUMNP)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

         CALL MAKCOL (LINK, NUMNP,
     &      NUMCOL, NEROW, A(KELCEN), A(KICOL1),
     &      A(KIRBOT), A(KIRTOP), A(KIBDIF),
     &      A(KNPBR), A(KNPBC), A(KNPTR), A(KNPTC))

         CALL MDDEL ('NPBROW')
         CALL MDDEL ('NPBCOL')
         CALL MDDEL ('NPTROW')
         CALL MDDEL ('NPTCOL')
         CALL MDDEL ('IRBDIF')

         CALL MDDEL ('ICOL1')

C      --Fill in the number of repetitions for the elements and nodes
C      --(NREL and NRNP) and the element column numbers (IELCOL) and
C      --the node rows and columns (NPCEN)

         CALL MDRSRV ('NRBOT', KNRBOT, NEROW)
         CALL MDRSRV ('NRTOP', KNRTOP, NEROW)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90

         N4 = NINT (DIM3 / 90)

         CALL FNPCEN (LINK, N4, NUMCOL, NEROW, A(KELCEN),
     &      A(KIRBOT), A(KIRTOP), NUMROW, NPCEN, NREL, IELCOL, NRNP,
     &      A(KNRBOT), A(KNRTOP))

         CALL MDDEL ('NRBOT')
         CALL MDDEL ('NRTOP')
         CALL MDDEL ('IRBOT')
         CALL MDDEL ('IRTOP')
         CALL MDDEL ('IELCEN')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 90
      END IF

C   --Calculate the new element numbers and increments

      IEL = 1
      JEL = 1
      DO 70 IELB = 1, NELBLK

         IF ((BLKTYP(IELB) .EQ. 'T') .OR. (BLKTYP(IELB) .EQ. 'S')) THEN

C         --Blocks that become multiple blocks do a plate/slice at a time

            DO 50 I = 1, NUMELB(IELB)
               IXEL(IEL) = JEL
               INCEL(IEL) = NUMELB(IELB)
               JEL = JEL + 1
               IEL = IEL + 1
   50       CONTINUE
            JEL = JEL + (NEREPL-1) * NUMELB(IELB)
         ELSE

C         --Other elements do an element at a time

            DO 60 I = 1, NUMELB(IELB)
               IXEL(IEL) = JEL
               INCEL(IEL) = 1
               JEL = JEL + NREL(IEL)
               IEL = IEL + 1
   60       CONTINUE
         END IF
   70 CONTINUE

C   --New number of elements
      NUMEL3 = JEL - 1

C   --Calculate the new node numbers

      JNP = 1
      DO 80 INP = 1, NUMNP
         IXNP(INP) = JNP
         JNP = JNP + NRNP(INP)
   80 CONTINUE

C   --New number of nodes
      NUMNP3 = JNP - 1

   90 CONTINUE
      RETURN
      END
