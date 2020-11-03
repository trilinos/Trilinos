C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE REMEL (IXELB, NUMLNK, LINK, VISELB,
     &                  IXNODE, IXELBO, IXELEM, NODIX, MAPEL)
C=======================================================================

      include 'ag_dbnums.blk'
      include 'ag_dbout.blk'
      include 'ag_remove.blk'

      INTEGER IXELB(0:NELBLK)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      LOGICAL VISELB(NELBLK)
      INTEGER IXNODE(*)
      INTEGER IXELBO(0:NELBLK)
      INTEGER IXELEM(*)
      INTEGER NODIX(*)
      INTEGER MAPEL(*)

      LOGICAL REMOVE

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

      IF (ISREMOVE) THEN

C      --Identify elements to be removed...

         NUMELO = 0
         IXLAST = IXELBO(0)
C        Loop 1 to number of element blocks
         DO IELB = 1, NELBLK
           DO IX = IXLAST+1, IXELBO(IELB)
C              iel = index of element
             IEL = IXELEM(IX)
             if (idsglobal) then
               IGEL = MAPEL(IEL)
             else
               igel = iel
             end if
             remove = .false.
             do ir = 1, irmcnt
               if (igel .eq. idsrem(ir)) then
                 write (*,999) iel, mapel(iel)
 999             FORMAT(' Removed element with local id ',I10,
     *             ', global id ',I10)
                 remove = .true.
                 go to 99
               end if
             end do
 99          continue
             if (.not. remove) then
               numelo = numelo + 1
               ixelem(numelo) = iel
             end if
           end do
           IXLAST = IXELBO(IELB)
           IXELBO(IELB) = NUMELO
         end do
      END IF

C   --Set the selected element blocks in the zoom mesh

      IF (ISREMOVE) THEN
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

