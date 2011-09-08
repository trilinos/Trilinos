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
      SUBROUTINE CKEB1 (IEL0, IELB, IDELB, NUMELB, NUMLNK, NUMNP, LINK,
     *  NODUSE)
C=======================================================================

C   --*** CKEB1 *** (GROPE) Check database element block connectivity
C   --
C   --CKEB1 checks that the database element block connectivity is within
C   --the nodal range.
C   --
C   --Parameters:
C   --   IEL0 - IN - the number of the first element in this block - 1
C   --   IELB - IN - the number of this element block
C   --   IDELB - IN - the element block ID for this block
C   --   NUMELB - IN - the number of elements for this block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMNP - IN - the number of nodes
C   --   LINK - IN - the connectivity array for this block
C   --   NODUSE - OUT - scratch array to determine whether all nodes used by an element.
      
      include 'errcnt.blk'
      INTEGER LINK(NUMLNK,*)
      INTEGER NODUSE(*)
      CHARACTER*132 STRA

      IF (NUMELB .LE. 0) GOTO 110
      IF (NUMLNK .LE. 0) GOTO 110

      NERR = 0

      DO 100 NE = 1, NUMELB
         CALL CHKRNG (LINK(1,NE), NUMLNK, NUMNP, NZERO, IERR)
         IF (IERR .GT. 0 .OR. NZERO .GT. 0) THEN
            IF (NERR .EQ. 0) THEN
               WRITE (*, 10000, IOSTAT=IDUM)
     &            'Connectivity Problems', IELB, IDELB
            END IF
            if (nerr .lt. maxerrs .or. maxerrs .le. 0) then
               WRITE (*, 10010, IOSTAT=IDUM)
     &              NE, NE+IEL0, (LINK(I,NE), I=1,NUMLNK)
            else if (nerr .eq. maxerrs .and. maxerrs .gt. 0) then
               call prterr('CMDSPEC',
     $              '...skipping additional errors...')
            end if
            NERR = NERR + 1
         END IF

         DO 90 I=1,NUMLNK
           NODE = LINK(I,NE)
           IF (NODE .LE. NUMNP) THEN
C ... Note that if NODE is out of range, the error message above will have
C     already been printed, so we dont print anything here.
             NODUSE(NODE) = 1
           END IF
 90      CONTINUE

 100  CONTINUE
      if (nerr .gt. 0) then
         write (stra, 10020) nerr, idelb
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
      end if
      

  110 CONTINUE
      RETURN

10000  FORMAT (/, 1X, '     #  elem      ', A, ' for block #', I10,
     &   ',', I10, ' = ID')
10010  FORMAT (1X, I10, I10, 5X, 6I10, :, /,
     &   (26X, 6I10))
10020    FORMAT('ELEMENT CONNECTIVITY ERROR: Found ',I10,
     $      ' errors in element connectivity check for element block '
     $      , i10)
      END
