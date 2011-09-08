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
      SUBROUTINE PRESS (OPTION, NOUT, NUMESS, LISESS, LESSEL, LESSNL,
     &  IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTSESS, FACESS,
     *  NAME,  nvar, namvar, isvok, lisvar, MAPEL, DOMAPE)
C=======================================================================

C     --*** PRESS *** (GROPE) Display database element side set
C     --
C     --PRESS displays the element side sets.
C     --
C     --Parameters:
C     --   OPTION - IN - '*' to print all, else print options:
C     --      ' ' - set summary (number of nodes, etc)
C     --      'E' - elements in set (may not be combined with 'N' or 'F')
C     --      'N' - nodes in set (may not be combined with 'E')
C     --      'F' - distribution factors for set (may not be combined with 'E')
C     --   NOUT - IN - the output file, <=0 for standard
C     --   NUMESS - IN - the number of element side sets
C     --   LISESS - IN - the indices of the selected element side sets
C     --   LESSEL - IN - the number of elements for all sets
C     --   LESSNL - IN - the number of nodes for all sets
C     --   IDESS - IN - the element side set ID for each set
C     --   NEESS - IN - the number of elements for each set
C     --   NNESS - IN - the number of nodes for each set
C     --   IXEESS - IN - the index of the first element for each set
C     --   IXNESS - IN - the index of the first node for each set
C     --   LTEESS - IN - the elements for all sets
C     --   LTSESS - IN - the element sides for all sets
C     --   FACESS - IN - the distribution factors for all sets
C     --   NVAR - IN - the number of variables
C     --   NAMVAR - IN - the names of the variables
C     --   ISVOK  - IN - the variable truth table;
C     --      variable i of set j exists iff ISVOK(i,j) is NOT 0
C     --   LISVAR  - SCRATCH - size = NVAR (if 'V' in OPTION)

      include 'params.blk'
      CHARACTER*(*) OPTION
      INTEGER LISESS(0:*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL FACESS(*)
      CHARACTER*(*) NAME(*)
      CHARACTER*(*) NAMVAR(*)
      INTEGER ISVOK(NVAR,*)
      INTEGER LISVAR(*)
      INTEGER MAPEL(*)
      LOGICAL DOMAPE

      LOGICAL ALLSAM
      LOGICAL DOELE, DONOD, DOFAC, DOVTBL
      CHARACTER*20 STRA, STRB, STRC

      DOELE  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0))
      DONOD  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOFAC  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))

      IF (NOUT .GT. 0) THEN
        IF (DOELE) THEN
          WRITE (NOUT, 10020) 'ELEMENT LIST'
        ELSE IF (DONOD .AND. DOFAC) THEN
          WRITE (NOUT, 10020) 'ELEM.FACE LIST AND DISTRIBUTION FACTORS'
        ELSE IF (DONOD) THEN
          WRITE (NOUT, 10020) 'ELEM.FACE LIST'
        ELSE IF (DOFAC) THEN
          WRITE (NOUT, 10020) 'DISTRIBUTION FACTORS'
        ELSE
          WRITE (NOUT, 10020)
        END IF
        if (domape) then
          write (nout, 10025)
        end if
      ELSE
        if (domape) then
          write (*, 10025)
        end if
      END IF

      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, *)
      ELSE
        WRITE (*, *)
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NUMESS
10000 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LESSEL
10010 FORMAT ('(index=', I10, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)
      WRITE (STRC, 10010, IOSTAT=IDUM) LESSNL
      CALL PCKSTR (1, STRC)
      LSTRC = LENSTR (STRC)

      DO 100 IX = 1, LISESS(0)
        IESS = LISESS(IX)
        WRITE (STRA, 10000, IOSTAT=IDUM) IESS
        CALL PCKSTR (1, STRA)
        WRITE (STRB, 10010, IOSTAT=IDUM) IXEESS(IESS)
        CALL PCKSTR (1, STRB)
        WRITE (STRC, 10010, IOSTAT=IDUM) IXNESS(IESS)
        CALL PCKSTR (1, STRC)
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10030, IOSTAT=IDUM)
     &      IDESS(IESS), STRA(:LSTRA),
     &      NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $      NAME(IESS)(:LENSTR(NAME(IESS)))
        ELSE
          WRITE (*, 10030, IOSTAT=IDUM)
     &      IDESS(IESS), STRA(:LSTRA),
     &      NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $      NAME(IESS)(:LENSTR(NAME(IESS)))
        END IF

        IF (DOELE .AND. (NEESS(IESS) .GT. 0)) THEN
          IS = IXEESS(IESS)
          IE = IS + NEESS(IESS) - 1
          IF (NOUT .GT. 0) THEN
            if (domape) then
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            else 
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (LTEESS(I), I=IS,IE)
            end if
          ELSE
            if (domape) then
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            else
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (LTEESS(I), I=IS,IE)
            end if
          END IF
        END IF

C     ... This used to print the nodes of the sideset faces, it now prints
C     the local element faces of the sideset         
        IF (DONOD .AND. (NEESS(IESS) .GT. 0)) THEN
          IS = IXEESS(IESS)
          IE = IS + NEESS(IESS) - 1
          IF (NOUT .GT. 0) THEN
            if (domape) then
              WRITE (NOUT, 10045, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)),LTSESS(I), I=IS,IE)
            else
              WRITE (NOUT, 10045, IOSTAT=IDUM)
     &          (LTEESS(I),LTSESS(I), I=IS,IE)
            endif
          ELSE
            if (domape) then
              WRITE (*, 10045, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)),LTSESS(I), I=IS,IE)
            else
              WRITE (*, 10045, IOSTAT=IDUM)
     &          (LTEESS(I),LTSESS(I), I=IS,IE)
            end if
          END IF
        END IF

        IF (DOVTBL) THEN
          NSEL = 0
          mxnam = 0
          DO 30 I = 1, NVAR
            IF (ISVOK(I,IESS) .NE. 0) THEN
              NSEL = NSEL + 1
              LISVAR(NSEL) = I
              LNAM = lenstr(namvar(i))
              mxnam = max(lnam, mxnam)
            END IF
 30       CONTINUE
          
          if (nsel .gt. 0) then
            IF (NOUT .GT. 0) THEN
              if (mxnam .gt. 24) then
                WRITE (NOUT, 10090, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              else 
                WRITE (NOUT, 10070, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              endif
              write (nout, 10080)
            ELSE
              if (mxnam .gt. 24) then
                WRITE (*, 10090, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              else 
                WRITE (*, 10070, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              endif
              write (*, 10080)
            END IF
          end if
        END IF
        
        IF (DOFAC) THEN
          IF (NNESS(IESS) .GT. 0) THEN
            IS = IXNESS(IESS)
            IE = IS + NNESS(IESS) - 1
C     ... See if all values are the same
            val = facess(is)
            allsam = .TRUE.
            do 50 i=is+1, ie
              if (facess(i) .ne. val) then
                allsam = .FALSE.
                go to 90
              end if
 50         continue
 90         continue
            if (allsam) then
              IF (NOUT .GT. 0) THEN
                WRITE (NOUT, 10055, IOSTAT=IDUM) VAL
              ELSE
                WRITE (*, 10055, IOSTAT=IDUM) VAL
              END IF
            else
C ... If it looks like the sideset elements are homogenous, then print
C     element/face/df information; otherwise, just print df...
              numdf = nness(iess)
              numel = neess(iess)
              ndfpe = numdf/numel
              if (ndfpe*numel .eq. numdf .and. numdf .gt. 0) then
                ISE = IXEESS(IESS)
                IEE = ISE + NEESS(IESS) - 1
                IDS = IS
                do i=ise, iee
                  if (domape) then
                    iel = mapel(lteess(i))
                  else
                    iel = lteess(i)
                  end if
                  if (nout .gt. 0) then
                    write (nout, 10100) iel, ltsess(i),
     *                (facess(j),j=ids,ids+ndfpe-1)
                  else
                    write (*, 10100) iel, ltsess(i),
     *                (facess(j),j=ids,ids+ndfpe-1)
                  end if
                  ids = ids + ndfpe
                end do
              else
                IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10050, IOSTAT=IDUM)
     &              (FACESS(I), I=IS,IE)
                ELSE
                  WRITE (*, 10050, IOSTAT=IDUM)
     &              (FACESS(I), I=IS,IE)
                END IF
              end if
            end if
          else
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10060, IOSTAT=IDUM)
            ELSE
              WRITE (*, 10060, IOSTAT=IDUM)
            END IF
          end if
        END IF
 100  CONTINUE
      
      RETURN

10020 FORMAT (/, 1X, 'ELEMENT SIDE SETS', :, ' - ', A)
10025 FORMAT (1x, 'Element Ids are Global')
10030 FORMAT (1X, 'Set', I10, 1X, A, ':',
     &  I10, ' elements', 1X, A,
     &  I10, ' nodes', 1X, A, ' name = "',A,'"')
10040 FORMAT ((1X, 8I10))
10045 FORMAT ((1X, 8(I10,'.',I1)))
10050 FORMAT ((1X, 6 (1X, 1pE11.4)))
10055 FORMAT (10x, 'All distribution factors are equal to ', 1pe11.4)
10060 FORMAT (10x, 'Distribution factors not stored in file.')
10070 FORMAT ((2x,4(2X, A)))
10090 FORMAT ((2x,3(2X, A)))
10080 FORMAT (1X)
10100 FORMAT ((1X, I10,'.',I1,8x,(8(1x,1pe11.4))))
      END
