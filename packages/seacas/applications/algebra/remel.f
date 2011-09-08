C    Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
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
      SUBROUTINE REMEL (IXELB, NUMLNK, LINK, VISELB,
     &                  IXNODE, IXELBO, IXELEM, NODIX, MAPEL)
C=======================================================================

      include 'dbnums.blk'
      include 'dbout.blk'
      include 'remove.blk'

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
C                 Indicies of the zoomed elements
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

