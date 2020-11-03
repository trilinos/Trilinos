C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FILTEL (IXELB, NUMLNK, LINK, VISELB,
     &                   IXNODE, IXELBO, IXELEM, NODIX,
     *                   ISEVOK, VALUES)
C=======================================================================

      include 'ag_dbnums.blk'
      include 'ag_dbout.blk'
      include 'ag_filter.blk'

      INTEGER IXELB(0:NELBLK)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      LOGICAL VISELB(NELBLK)
      INTEGER IXNODE(*)
      INTEGER IXELBO(0:NELBLK)
      INTEGER IXELEM(*)
      INTEGER NODIX(*)
      LOGICAL ISEVOK(*)
      REAL    VALUES(*)

      DO 100 INP = 1, NUMNP
        NODIX(INP) = 0
 100  CONTINUE

C     If (Output Element Blocks < Element Blocks)
      IF (NUMEQL (.TRUE., NELBLK, VISELB) .LT. NELBLK) THEN

C      --Identify elements of a selected element block, and save index
C      --Eliminate nodes that are not in an element of a selected element block

         NUMELO = 0
         IXELBO(0) = NUMELO
C        Loop 1 to number of element blocks
         DO 130 IELB = 1, NELBLK
C           Write element block IELB
            IF (VISELB(IELB)) THEN
C              0 index into cumulative count array
C              Loop over cumulative number of elements in element block IELB
               DO 120 IEL = IXELB(IELB-1)+1, IXELB(IELB)
C                 increment number of elements to output
                  NUMELO = NUMELO + 1
C                 Indices of the zoomed elements
                  IXELEM(NUMELO) = IEL
  120          CONTINUE
            END IF
C           Set cumul element count for output block IELB
            IXELBO(IELB) = NUMELO
  130    CONTINUE

      ELSE
         NUMELO = 0
         IXELBO(0) = NUMELO
C        loop from 1 to number of element blocks
         DO 160 IELB = 1, NELBLK
C           Loop over cumulative elements count in element block IELB
            DO 150 IEL = IXELB(IELB-1)+1, IXELB(IELB)
C              Increment number of elements to output
               NUMELO = NUMELO + 1
C              Initialize indices of the zoomed elements
               IXELEM(NUMELO) = IEL
  150       CONTINUE
C           Set the number of elements to output for element block IELB
            IXELBO(IELB) = NUMELO
  160    CONTINUE
      END IF

      IF (ISFILTER) THEN

C      --Identify elements to be retained...

         NUMELO = 0
         IXLAST = IXELBO(0)
C        Loop 1 to number of element blocks
         DO IELB = 1, NELBLK
C           Loop over cumulative element count in element block IELB
           if (isevok(ielb) .and. viselb(ielb)) then
             if (cmpflt .eq. 1) then
               call dolt(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             else if (cmpflt .eq. 2) then
               call dole(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             else if (cmpflt .eq. 3) then
               call doeq(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             else if (cmpflt .eq. 4) then
               call done(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             else if (cmpflt .eq. 5) then
               call dogt(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             else if (cmpflt .eq. 6) then
               call doge(ixlast+1, ixelbo(ielb), ixelem, values,
     *           valflt, numelo)
             end if
           else
C ... If the variable doesn't exist on a block, retain all elements in the block
            DO IX = IXLAST+1, IXELBO(IELB)
C              iel = index of element
               IEL = IXELEM(IX)
               numelo = numelo + 1
               ixelem(numelo) = iel
             end do
           end if
           IXLAST = IXELBO(IELB)
           IXELBO(IELB) = NUMELO
         end do
      END IF

C   --Set the selected element blocks in the zoom mesh

      IF (ISFILTER) THEN
         DO 220 IELB = 1, NELBLK
            NELO = IXELBO(IELB) - IXELBO(IELB-1)
            VISELB(IELB) = NELO .GT. 0
  220    CONTINUE
      END IF

C   --Count the number of selected element blocks

      NELBO = NUMEQL (.TRUE., NELBLK, VISELB)

C   --Index nodes of selected element block within the zoom mesh
      NUMELO = 0
      IXLAST = IXELBO(0)
C        Loop 1 to number of element blocks
      DO IELB = 1, NELBLK
C       Loop over cumulative element count in element block IELB
        DO IX = IXLAST+1, IXELBO(IELB)
C         iel = index of zoomed element
          IEL = IXELEM(IX)
          numelo = numelo + 1
C         ixlo = 0 index for elemtn block
          IXL0 = IDBLNK (IELB, IEL, IXELB, NUMLNK) - 1
C         Loop 1 to number of nodes per element - element block IELB
          DO K = 1, NUMLNK(IELB)
            NODIX(LINK(IXL0+K)) = -1
          end do
        end do
        IXLAST = IXELBO(IELB)
        IXELBO(IELB) = NUMELO
      end do

      NUMNPO = 0
      DO INP = 1, NUMNP
        IF (NODIX(INP) .NE. 0) THEN
          NUMNPO = NUMNPO + 1
          IXNODE(NUMNPO) = INP
          NODIX(INP) = NUMNPO
        END IF
      end do

      WRITE (*, 10000) NUMNPO, NUMELO, NELBO
10000  FORMAT (
     &   /, 1X, 'Number of output nodes               =', I10
     &   /, 1X, 'Number of output elements            =', I10
     &   /, 1X, 'Number of output element blocks      =', I10
     &   )

      RETURN
      END

      subroutine dolt(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .lt. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end

      subroutine dole(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .le. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end

      subroutine doeq(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .eq. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end

      subroutine done(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .ne. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end

      subroutine dogt(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .gt. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end

      subroutine doge(ibeg, iend, ixelem, values, valflt, numelo)
      integer ixelem(*)
      real values(*)

      DO IX = ibeg, iend
C              iel = index of element
        IEL = IXELEM(IX)
        if (.not.(values(iel) .ge. valflt)) then
C ... retain element
          numelo = numelo + 1
          ixelem(numelo) = iel
        end if
      end do
      return
      end
