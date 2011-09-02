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
      SUBROUTINE PREB1 (OPTION, NOUT, IEL0, NLISEL, LISEL,
     &  NLINK, NATR, LINK, ATRIB, natrdm, 
     *  MAPN, DOMAPN, MAPE, DOMAPE)
C=======================================================================

C     --*** PREB1 *** (BLOT) Display database element block information
C     --
C     --PREB1 displays the element block internal information (either the
C     --connectivity or the attributes.
C     --
C     --Parameters:
C     --   OPTION - IN - '*' to print all, else print options:
C     --      'C' to print connectivity
C     --      'A' to print attributes
C     --   NOUT - IN - the output file, <=0 for standard
C     --   IEL0 - IN - the element offset of the elements in the block
C     --   NLISEL - IN - the number of selected elements
C     --   LISEL - IN - the indices of the selected elements
C     --   NLINK - IN - the number of nodes per element
C     --   NATR - IN - the number of attributes per element
C     --   LINK - IN - the connectivity array for this block
C     --   ATRIB - IN - the attribute array for this block

      CHARACTER*(*) OPTION
      INTEGER LISEL(*)
      INTEGER LINK(NLINK,*)
      REAL ATRIB(natrdm,*)
      INTEGER MAPN(*), MAPE(*)
      LOGICAL DOMAPN, DOMAPE

      LOGICAL ALLSAM
      LOGICAL DOCONN, DOATR

      IF (NLISEL .LE. 0) RETURN

      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
     &  .AND. (NLINK .GT. 0)
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))
     &  .AND. (NATR .GT. 0)
      IF (.NOT. (DOCONN .OR. DOATR)) RETURN

      if (doconn) then
        if (nout .gt. 0) then
          WRITE (NOUT, 10000) 'Connectivity'
        else
          WRITE (*, 10000) 'Connectivity'
        end if
        do 100 ix=1, nlisel
          IEL = LISEL(IX)
          if (domape) then
            ID = MAPE(IEL)
          else
            id = iel
          end if
          NE = IEL - IEL0
          IF (NOUT .GT. 0) THEN
            if (domapn) then
              WRITE (NOUT, 10010, IOSTAT=IDUM)
     &          NE, ID, (MAPN(LINK(I,NE)), I=1,NLINK)
            else
              WRITE (NOUT, 10010, IOSTAT=IDUM)
     &          NE, ID, (LINK(I,NE), I=1,NLINK)
            end if
          ELSE
            if (domapn) then
              WRITE (*, 10010, IOSTAT=IDUM)
     &          NE, ID, (MAPN(LINK(I,NE)), I=1,NLINK)
            else
              WRITE (*, 10010, IOSTAT=IDUM)
     &          NE, ID, (LINK(I,NE), I=1,NLINK)
            end if
          END IF
 100    continue
      END IF
      
      if (doatr) then
C     ... See if all attributes are the same
        allsam = .true.
        do 180 ia = 1, natr
          iel = lisel(1)-iel0
          aval = atrib(ia, iel)
          do 170 ix = 2, nlisel
            iel = lisel(ix) - iel0
            if (aval .ne. atrib(ia, iel)) then
              allsam = .false.
              go to 190
            end if
 170      continue
 180    continue

C     ... Print either all values (if not all same), or the common values
 190    continue
        if (allsam) then
          ne = lisel(1) - iel0
          IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10025, IOSTAT=IDUM)
     &        (ATRIB(I,NE), I=1,NATR)
          ELSE
            WRITE (*, 10025, IOSTAT=IDUM)
     &        (ATRIB(I,NE), I=1,NATR)
          END IF
        else
          if (nout .gt. 0) then
            WRITE (NOUT, 10000) 'Attributes'
          else
            WRITE (*, 10000) 'Attributes'
          end if
          DO 200 IX = 1, NLISEL
            IEL = LISEL(IX)
            NE = IEL - IEL0
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10020, IOSTAT=IDUM)
     &          NE, IEL, (ATRIB(I,NE), I=1,NATR)
            ELSE
              WRITE (*, 10020, IOSTAT=IDUM)
     &          NE, IEL, (ATRIB(I,NE), I=1,NATR)
            END IF
 200      CONTINUE
        end if
      end if
      
      RETURN

10000 FORMAT (1X, '         #      elem            ', A)
10010 FORMAT (1X, I10, I10, 5X, 8I10, :, /,
     &  (26X, 8I10))
10020 FORMAT (1X, I10, I10, 3X, 4 (1X, 1pE11.4), :, /,
     &  (24X, 4 (1X, 1pE11.4)))
10025 FORMAT (2X, 'All attributes are equal. Values are: ',
     &  (3 (1x, 1pE11.4)), :, /,
     &  (3X, 6 (1X, 1pE11.4)))
10030 FORMAT (1x, 'Ids are global')
      END
