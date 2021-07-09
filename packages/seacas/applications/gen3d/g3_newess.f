C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWESS (IDFRO, IDBCK, LINK,
     &   ISSFRO, ISSBCK, NSSUR, NESUR, NSSFRO, NSSBCK,
     &   IDESS, NEESS, NEES3, NNESS, NNES3,
     &   IXEESS, IXEES3, IXNESS, IXNES3,
     &   LTEESS, LTEES3, LTSESS, LTSES3, LTNESS, LTNES3, FACESS, FACES3,
     &   IXEL, INCEL, NREL, IELCOL, IXNP, NRNP, IROT)
C=======================================================================

C   --*** NEWESS *** (GEN3D) Calculate 3D side sets
C   --   Written by Amy Gilkey - revised 01/12/88
C   --
C   --NEWESS calculates the side set information for the 3D database.
C   --The elements in the front and back side sets are not calculated
C   --(since they are easily derived), and they are not included in the
C   --length tally.  The nodes in the front and back node sets are
C   --calculated, and their number is returned, but they are not included
C   --in the length tally.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface side sets; (0) = length
C   --   IDBCK - IN - ids for back surface side sets; (0) = length
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   ISSFRO - OUT - the elements in the front surface side set
C   --   ISSBCK - OUT - the elements in the back surface side set
C   --   NSSUR - OUT - the number of nodes in the surface side set
C   --   NESUR - OUT - the number of elements in the surface side set
C   --   NSSFRO - OUT - the nodes in the front surface side set
C   --   NSSBCK - OUT - the nodes in the back surface side set
C   --   IDESS - IN - the ids for each 2D set
C   --   NEESS - IN - the number of elements for each 2D set
C   --   NEES3 - OUT - the number of elements for each 3D set
C   --   NNESS - IN - the number of nodes for each 2D set
C   --   NNES3 - OUT - the number of nodes for each 3D set
C   --   IXEESS - IN - the index of the first element for each 2D set
C   --   IXEES3 - OUT - the index of the first element for each 3D set
C   --   IXNESS - IN - the index of the first node for each 2D set
C   --   IXNES3 - OUT - the index of the first node for each 3D set
C   --   LTEESS - IN - the elements for all 2D sets
C   --   LTEES3 - OUT - the elements for all 3D sets
C   --   LTSESS - IN - the element sides for all 2D sets
C   --   LTSES3 - OUT - the element sides for all 3D sets
C   --   LTNESS - IN - the nodes for all 2D sets
C   --   LTNES3 - OUT - the nodes for all 3D sets
C   --   FACESS - IN - the distribution factors for all 2D sets
C   --   FACES3 - OUT - the distribution factors for all 3D sets
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   IROT - IN - offset from original link order
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Sets LESSEO, LESSNO of /DBNUM3/
C   --   Uses NNREPL, NEREPL, DIM3 of /PARAMS/

      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER LINK(4,*)
      INTEGER ISSFRO(NUMEL), ISSBCK(NUMEL)
      INTEGER NSSFRO(*), NSSBCK(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*), NEES3(*)
      INTEGER NNESS(*), NNES3(*)
      INTEGER IXEESS(*), IXEES3(*)
      INTEGER IXNESS(*), IXNES3(*)
      INTEGER LTEESS(*), LTEES3(*)
      INTEGER LTSESS(*), LTSES3(*)
      INTEGER LTNESS(*), LTNES3(*)
      REAL FACESS(*), FACES3(*)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)
      INTEGER IXNP(*), NRNP(*), IROT(*)

      LOGICAL ROT360
      LOGICAL SWAPIT

      ROT360 = (NEREPL .EQ. NNREPL)
      N4 = NINT (DIM3 / 90)

C   --Side set ids - unchanged

      CONTINUE

C   --Side set elements - 2D array must be different from 3D array to prevent
C   --overwriting; add on elements for each plate/slice

C   --Side set nodes and distribution factors - 2D array must be different
C   --from 3D array to prevent overwriting; first swap 2D side set nodes;
C   --then complete 2D side set by adding nodes from next plate/slice in
C   --reverse order (2-1-5-6), then add on nodes for each plate/slice

C   --Swap nodes to reverse order

      DO 30 IESS = 1, NUMESS
         IP0 = IXNESS(IESS) - 1
         NES = 2
         DO 20 J = 1, NEESS(IESS)
            IDOWN = NES
            DO 10 I = 1, NES/2
               ISV = LTNESS(IP0+I)
               LTNESS(IP0+I) = LTNESS(IP0+IDOWN)
               LTNESS(IP0+IDOWN) = ISV
               IDOWN = IDOWN - 1
   10       CONTINUE
            IP0 = IP0 + NES
   20    CONTINUE
   30 CONTINUE

      NE = 0
      N = 0
      DO 130 IESS = 1, NUMESS
C      --Index of each side set elements
         IXEES3(IESS) = NE + 1
C      --Index of each side set nodes
         IXNES3(IESS) = N + 1

         IP0 = IXNESS(IESS) - 1
         IXE0 = IXEESS(IESS) - 1
         NES = 2

         DO 120 J = 1, NEESS(IESS)
            IEL = LTEESS(IXE0+J)
C ... If IROT(IEL) .ne. 1, then the connectivity of the 2D element
C     was rotated in ixlink.  Update the side number to match the
C     new connectivity.
            ISD = LTSESS(IXE0+J) - IROT(IEL) + 1
            if (isd .le. 0) isd = isd + 4
            JEL = IXEL(IEL)

            IF (IELCOL(IEL) .EQ. 0) THEN
               DO 40 NR = 1, NREL(IEL)
                  NE = NE + 1
                  LTEES3(NE) = JEL
                  LTSES3(NE) = ISD
                  JEL = JEL + INCEL(IEL)
   40          CONTINUE

               DO 60 NR = 1, NREL(IEL)
                  NEND = N + 2*NES + 1
                  DO 50 I = 1, NES
                     INP = LTNESS(IP0+I)
                     JNP = IXNP(INP)
                     N = N + 1
                     LTNES3(N) = JNP + NR-1
                     IF (NR .LT. NNREPL) THEN
                        LTNES3(NEND-I) = JNP + NR
                     ELSE
                        LTNES3(NEND-I) = JNP
                     END IF
                     FACES3(N) = FACESS(IP0+I)
                     FACES3(NEND-I) = FACESS(IP0+I)
   50             CONTINUE
                  N = N + NES
   60          CONTINUE

            ELSE

C            --Handle center element

               NE4 = NREL(IEL) / N4
               NCORN = INT (NE4/2) + 1

               INP1 = LTNESS(IP0+1)
               INP2 = LTNESS(IP0+2)

               IF (NRNP(INP1) .EQ. NRNP(INP2)) THEN

                  DO 70 NR = 1, NREL(IEL)
                     IF (NR .EQ. NCORN) THEN
                        IF (NRNP(INP2) .GT. NREL(IEL)) THEN
                           NE = NE + 1
                           LTEES3(NE) = JEL
                           LTSES3(NE) = ISD
                           NE = NE + 1
                           LTEES3(NE) = JEL
                           LTSES3(NE) = ISD
                        ELSE
                           CONTINUE
                        END IF
                        NCORN = NCORN + NE4
                     ELSE
                        NE = NE + 1
                        LTEES3(NE) = JEL
                        LTSES3(NE) = ISD
                     END IF
                     JEL = JEL + INCEL(IEL)
   70             CONTINUE

                  DO 90 NR = 1, NRNP(INP2)
                     NEND = N + 2*NES + 1
                     DO 80 I = 1, NES
                        INP = LTNESS(IP0+I)
                        JNP = IXNP(INP)
                        N = N + 1
                        LTNES3(N) = JNP + NR-1
                        FACES3(N) = FACESS(IP0+I)
                        IF (.NOT. ROT360 .OR. NR .LT. NRNP(INP2)) THEN
                           LTNES3(NEND-I) = JNP + NR
                        ELSE
                           LTNES3(NEND-I) = JNP
                        END IF
                        FACES3(NEND-I) = FACESS(IP0+I)
   80                CONTINUE
                     N = N + NES
   90             CONTINUE

               ELSE
                  DO 100 NR = 1, NREL(IEL)
                     NE = NE + 1
                     LTEES3(NE) = JEL
                     LTSES3(NE) = ISD
                     JEL = JEL + INCEL(IEL)
  100             CONTINUE

                  SWAPIT = (NRNP(INP1) .GT. NRNP(INP2))
                  IF (SWAPIT) THEN
                     I = INP1
                     INP1 = INP2
                     INP2 = I
                  END IF
                  JNP1 = IXNP(INP1)
                  JNP2 = IXNP(INP2)
                  I1 = 0
                  I2 = 0
                  DO 110 NR = 1, NREL(IEL)
                     LTNES3(N+1) = JNP1 + I1
                     LTNES3(N+2) = JNP2 + I2
                     IF (NR .EQ. NCORN) THEN
                        I3 = I2 + 2
                        IF (ROT360 .AND. I3 .GE. NRNP(INP2)) I3 = 0
                        LTNES3(N+4) = JNP2 + I3
                     ELSE
                        I3 = I1 + 1
                        IF (ROT360 .AND. I3 .GE. NRNP(INP1)) I3 = 0
                        LTNES3(N+4) = JNP1 + I3
                     END IF
                     I4 = I2 + 1
                     IF (ROT360 .AND. I4 .GE. NRNP(INP2)) I4 = 0
                     LTNES3(N+3) = JNP2 + I4
                     IF (SWAPIT) THEN
                        I = LTNES3(N+1)
                        LTNES3(N+1) = LTNES3(N+2)
                        LTNES3(N+2) = I
                        I = LTNES3(N+3)
                        LTNES3(N+3) = LTNES3(N+4)
                        LTNES3(N+4) = I
                     END IF
                     FACES3(N+1) = FACESS(IP0+1)
                     FACES3(N+2) = FACESS(IP0+2)
                     FACES3(N+3) = FACESS(IP0+2)
                     FACES3(N+4) = FACESS(IP0+1)
                     IF (NR .EQ. NCORN) THEN
                        NCORN = NCORN + NE4
                        I1 = I1
                        I2 = I2 + 2
                     ELSE
                        I1 = I1 + 1
                        I2 = I2 + 1
                     END IF
                     N = N + 4
  110             CONTINUE
               END IF
            END IF

            IP0 = IP0 + NES
  120    CONTINUE

C      --Number of elements in each set
         NEES3(IESS) = NE - IXEES3(IESS) + 1
C      --Number of nodes in each set
         NNES3(IESS) = N - IXNES3(IESS) + 1
  130 CONTINUE

C   --Number of elements in all sets
      LESSEO = NE
C   --Number of nodes in all sets
      LESSNO = N

C   --Set up elements and nodes for front and back side sets

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)

      NESUR = 0
      NSSUR = 0
      IF ((NFRO .GT. 0) .OR. (NBCK .GT. 0)) THEN
         N = 0
         DO 150 IEL = 1, NUMEL
C ... Check if BAR or QUAD.  All are stored with 4-node connectivity
C     BAR will have last two connectivity entries equal to -1
C     The connectivity is initialized to -1 in rdelb.

           IF (LINK(3,IEL) .NE. -1 .AND. LINK(4,IEL) .NE. -1) then
             NESUR = NESUR + 1
             JEL = IXEL(IEL)
             IF (NFRO .GT. 0) ISSFRO(IEL) = JEL
             IF (NBCK .GT. 0) THEN
               ISSBCK(IEL) = JEL + nrel(iel)*INCEL(IEL)-1
C ... The center element on the back side needs special treatment when
C     the sideset faces are written out (in wress).  Flag it with
C     a negative element number (and remember to change back in wress)
               if (iscent .and. nrel(iel) .eq. 2)
     *           issbck(iel) = -issbck(iel)
             end if
             NLINK = 4
             NEND = N + NLINK + 1
             DO I = 1, NLINK
               INP = LINK(I,IEL)
               JNP = IXNP(INP)
               N = N + 1
               IF (NFRO .GT. 0) NSSFRO(N) = JNP
               IF (NBCK .GT. 0) NSSBCK(NEND-I) = JNP + NRNP(INP)-1
             end do
           ELSE IF (LINK(3,IEL) .NE. -1 .AND. LINK(4,IEL) .EQ. -1) then
C ... Triangle to Wedge
             NESUR = NESUR + 1
             JEL = IXEL(IEL)
C ... Need to tell wress that this is a wedge and not a hex.
C     Use the same negative element kluge since the faces
C     work out the same... Flag it with
C     a negative element number (and remember to change back in wress)
             IF (NFRO .GT. 0) ISSFRO(IEL) = -JEL
             IF (NBCK .GT. 0) then
               ISSBCK(IEL) = -(JEL + nrel(iel)*INCEL(IEL)-1)
             end if
             NLINK = 3
             NEND = N + NLINK + 1
             DO I = 1, NLINK
               INP = LINK(I,IEL)
               JNP = IXNP(INP)
               N = N + 1
               IF (NFRO .GT. 0) NSSFRO(N) = JNP
               IF (NBCK .GT. 0) NSSBCK(NEND-I) = JNP + NRNP(INP)-1
             end do
           end if

 150     CONTINUE
         NSSUR = N
       END IF

      RETURN
      END
