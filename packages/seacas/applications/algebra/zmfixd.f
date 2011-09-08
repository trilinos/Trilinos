C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C=======================================================================
      SUBROUTINE ZMFIXD (XN, YN, ZN, IXELB, NUMLNK, LINK, VISELB,
     &  IXNODE, IXELBO, IXELEM, NODIX)
C=======================================================================

C   --*** ZMFIXD *** (ALGEBRA) Zoom database fixed arrays
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --ZMFIXD finds the nodes within the zoomed coordinates and the elements
C   --within the zoomed mesh (if any nodes are in).  The connectivity
C   --array is fixed to point to the renumbered nodes.  The number of
C   --element blocks is recalculated.
C   --
C   --Parameters:
C   --   XN, YN, ZN, - IN - the coordinates
C   --   IXELB   - IN - the cumulative element counts for each element block
C   --   NUMLNK  - IN - the number of nodes per element in each element block
C   --   LINK    - IN - the connectivity for all elements
C   --   VISELB(i) - IN/OUT - true iff element block i is to be written
C   --   IXNODE   - OUT - the indices of the zoomed nodes (iff ISZOOM)
C   --   IXELBO   - OUT - the cumulative element counts for each output block
C   --   IXELEM   - OUT - the indices of the zoomed elements (iff ISZOOM)
C   --   NODIX(i) - OUT - the zoom mesh index for each node (iff ISZOOM)
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK of /DBNUMS/
C   --   Sets NUMNPO, NUMELO, NELBO of /DBOUT/
C   --   Uses ISZOOM, ZMLIM of /ZOOM/
      
      include 'dbnums.blk'
      include 'dbout.blk'
      include 'zoom.blk'

      REAL XN(*), YN(*), ZN(*)
      INTEGER IXELB(0:NELBLK)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      LOGICAL VISELB(NELBLK)
      INTEGER IXNODE(*)
      INTEGER IXELBO(0:NELBLK)
      INTEGER IXELEM(*)
      INTEGER NODIX(*)

      LOGICAL NODEIN, NODEOUT

C     If (Output Element Blocks < Element Blocks)
      IF (NUMEQL (.TRUE., NELBLK, VISELB) .LT. NELBLK) THEN

C      --Identify elements of a selected element block, and save index
C      --Eliminate nodes that are not in an element of a selected element block

C        initialize zoom mesh index for each node
        DO 100 INP = 1, NUMNP
          NODIX(INP) = 0
 100    CONTINUE

        NUMELO = 0
        IXELBO(0) = NUMELO
C        Loop 1 to number of element blocks
        DO 130 IELB = 1, NELBLK
C           Write element block IELB 
          IF (VISELB(IELB)) THEN
C              0 index into cumulative count array
            IXL0 = IDBLNK (IELB, 0, IXELB, NUMLNK) - 1
C              Loop over cumulative number of elements in element block IELB
            DO 120 IEL = IXELB(IELB-1)+1, IXELB(IELB)
C                 increment number of elements to output
              NUMELO = NUMELO + 1
C                 Indicies of the zoomed elements
              IXELEM(NUMELO) = IEL
C                 Loop 1 to number of nodes per element in elem blk IELB
              DO 110 K = 1, NUMLNK(IELB)
C                    initialize zoom mesh index for each node
                NODIX(LINK(IXL0+K)) = 1
 110          CONTINUE
C                 increment IXLO by the number of nodes per element block
              IXL0 = IXL0 + NUMLNK(IELB)
 120        CONTINUE
          END IF
C           Set cumul element count for output block IELB
          IXELBO(IELB) = NUMELO
 130    CONTINUE

      ELSE
C        initialize the zoom mesh index from 1 to number of nodes
        DO 140 INP = 1, NUMNP
          NODIX(INP) = 1
 140    CONTINUE

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
 150      CONTINUE
C           Set the number of elements to output for element block IELB
          IXELBO(IELB) = NUMELO
 160    CONTINUE
      END IF

      IF (ISZOOM) THEN

C      --Eliminate nodes outside the zoom mesh

C        Loop from 1 to total number of nodes
        DO  INP = 1, NUMNP
          IF (NODIX(INP) .GE. 1) THEN
            IF (NDIM .EQ. 2) THEN
              NODEIN =
     &          (XN(INP) .GE. ZMLIM(1)) .AND.
     &          (XN(INP) .LE. ZMLIM(2)) .AND.
     &          (YN(INP) .GE. ZMLIM(3)) .AND.
     &          (YN(INP) .LE. ZMLIM(4))
            ELSE IF (NDIM .GE. 3) THEN
              NODEIN =
     &          (XN(INP) .GE. ZMLIM(1)) .AND.
     &          (XN(INP) .LE. ZMLIM(2)) .AND.
     &          (YN(INP) .GE. ZMLIM(3)) .AND.
     &          (YN(INP) .LE. ZMLIM(4)) .AND.
     &          (ZN(INP) .GE. ZMLIM(5)) .AND.
     &          (ZN(INP) .LE. ZMLIM(6))
            ELSE
              NODEIN = .TRUE.
            END IF
C              Resest NODIX variable
            if (zoomin) then
              IF (.NOT. NODEIN) NODIX(INP) = 0
            else
              IF (NODEIN) NODIX(INP) = 0
            end if
          END IF
        END DO

C      --Identify elements within the zoom mesh, and save index

        NUMELO = 0
        IXLAST = IXELBO(0)
C        Loop 1 to number of element blocks
        DO 210 IELB = 1, NELBLK
C           Loop over cumulative element count in element block IELB
          DO 200 IX = IXLAST+1, IXELBO(IELB)
C              iel = index of zoomed element
            IEL = IXELEM(IX)
C              ixlo = 0 index for elemtn block
            IXL0 = IDBLNK (IELB, IEL, IXELB, NUMLNK) - 1
            N = 0
C              Loop 1 to number of nodes per element - element block IELB
            DO 180 K = 1, NUMLNK(IELB)
C                 If (NODIX(node#)>0) increment number of nodes to output
C                    for element block IELB
              IF (NODIX(LINK(IXL0+K)) .GT. 0) N = N + 1
 180        CONTINUE
            IF (N .GT. 0) THEN
C                 if (output nodes/element < nodes/element)
              IF (N .LT. NUMLNK(IELB)) THEN
C                    Loop 1 to number of nodes per element - element block IELB
                DO 190 K = 1, NUMLNK(IELB)
C                       if node is not to be output
C                          set nodix = -1
                  IF (NODIX(LINK(IXL0+K)) .EQ. 0)
     &              NODIX(LINK(IXL0+K)) = -1
 190            CONTINUE
              END IF
C                 Increment number of elements to output
              NUMELO = NUMELO + 1
              IXELEM(NUMELO) = IEL
            END IF
 200      CONTINUE
          IXLAST = IXELBO(IELB)
          IXELBO(IELB) = NUMELO
 210    CONTINUE
      END IF

C   --Set the selected element blocks in the zoom mesh

      IF (ISZOOM) THEN
        DO 220 IELB = 1, NELBLK
          NELO = IXELBO(IELB) - IXELBO(IELB-1)
          VISELB(IELB) = NELO .GT. 0
 220    CONTINUE
      END IF

C   --Count the number of selected element blocks

      NELBO = NUMEQL (.TRUE., NELBLK, VISELB)

C   --Index nodes of selected element block within the zoom mesh

      NUMNPO = 0
      DO 230 INP = 1, NUMNP
        IF (NODIX(INP) .NE. 0) THEN
          NUMNPO = NUMNPO + 1
          IXNODE(NUMNPO) = INP
          NODIX(INP) = NUMNPO
        END IF
 230  CONTINUE

      WRITE (*, 10000, IOSTAT=IDUM) NUMNPO, NUMELO, NELBO
10000 FORMAT (
     &  /, 1X, 'Number of output nodes               =', I10
     &  /, 1X, 'Number of output elements            =', I10
     &  /, 1X, 'Number of output element blocks      =', I10
     &  )

      RETURN
      END
