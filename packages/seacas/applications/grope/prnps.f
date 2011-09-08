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
      SUBROUTINE PRNPS (OPTION, NOUT, NUMNPS, LISNPS, LNPSNL,
     &     IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS, LTNNPS, FACNPS, NAME,
     $     nvar, namvar, isvok, lisvar, 
     $     MAPNO, DOMAPN)
C=======================================================================

C   --*** PRNPS *** (GROPE) Display database nodal point set
C   --
C   --PRNPS displays the nodal point sets.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' - set summary (number of nodes, etc)
C   --      'N' - nodes in set
C   --      'F' - distribution factors for set
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LISNPS - IN - the indices of the selected nodal point sets
C   --   LNPSNL - IN - the number of nodes for all sets
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   NDNPS - IN - the number of distribution factors for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   IXDNPS - IN - the index of the first dist factor for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --   NVAR - IN - the number of variables
C   --   NAMVAR - IN - the names of the variables
C   --   ISVOK  - IN - the variable truth table;
C   --      variable i of set j exists iff ISVOK(i,j) is NOT 0
C   --   LISVAR  - SCRATCH - size = NVAR (if 'V' in OPTION)

      include 'params.blk'
      CHARACTER*(*) OPTION
      INTEGER LISNPS(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXNNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      CHARACTER*(*) NAME(*)
      CHARACTER*(*) NAMVAR(*)
      INTEGER ISVOK(NVAR,*)
      INTEGER LISVAR(*)
      INTEGER MAPNO(*)
      LOGICAL DOMAPN

      LOGICAL ALLSAM
      LOGICAL DONOD, DOFAC, DOVTBL
      CHARACTER*20 STRA, STRB

      DONOD  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOFAC  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))

      IF (NOUT .GT. 0) THEN
        IF (DONOD .AND. DOFAC) THEN
          WRITE (NOUT, 10020) 'NODE LIST AND DISTRIBUTION FACTORS'
        ELSE IF (DONOD) THEN
          WRITE (NOUT, 10020) 'NODE LIST'
        ELSE IF (DOFAC) THEN
          WRITE (NOUT, 10020) 'DISTRIBUTION FACTORS'
        ELSE
          WRITE (NOUT, 10020)
        END IF
      END IF

      if (domapn) then
        if (nout .gt. 0) then
          write (nout, 10005)
        else
          write (*, 10005)
        end if
      end if
      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, *)
      ELSE
        WRITE (*, *)
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NUMNPS
10000 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LNPSNL
10010 FORMAT ('(index=', I10, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      DO 100 IX = 1, LISNPS(0)
        INPS = LISNPS(IX)
        WRITE (STRA, 10000, IOSTAT=IDUM) INPS
        CALL PCKSTR (1, STRA)
        WRITE (STRB, 10010, IOSTAT=IDUM) IXNNPS(INPS)
        CALL PCKSTR (1, STRB)
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10030, IOSTAT=IDUM)
     &      IDNPS(INPS), STRA(:LSTRA),
     &      NNNPS(INPS), STRB(:LSTRB),
     $          NAME(INPS)(:LENSTR(NAME(INPS)))
        ELSE
          WRITE (*, 10030, IOSTAT=IDUM)
     &      IDNPS(INPS), STRA(:LSTRA),
     &      NNNPS(INPS), STRB(:LSTRB),
     $          NAME(INPS)(:LENSTR(NAME(INPS)))
        END IF

        IF (DONOD .AND. (NNNPS(INPS) .GT. 0)) THEN
          IS = IXNNPS(INPS)
          IE = IS + NNNPS(INPS) - 1
          IF (NOUT .GT. 0) THEN
            if (domapn) then
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPNO(LTNNPS(I)), I=IS,IE)
            else
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (LTNNPS(I), I=IS,IE)
            end if
          ELSE
            if (domapn) then
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPNO(LTNNPS(I)), I=IS,IE)
            else
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (LTNNPS(I), I=IS,IE)
              end if
          END IF
        END IF

         IF (DOVTBL) THEN
            NSEL = 0
            mxnam = 0
            DO 30 I = 1, NVAR
               IF (ISVOK(I,INPS) .NE. 0) THEN
                  NSEL = NSEL + 1
                  LISVAR(NSEL) = I
                  LNAM = lenstr(namvar(i))
                  mxnam = max(lnam, mxnam)
               END IF
 30         CONTINUE
            
            if (nsel .gt. 0) then
              IF (NOUT .GT. 0) THEN
                if (mxnam .gt. 24) then
                  WRITE (NOUT, 10090, IOSTAT=IDUM)
     &              (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
                else 
                  WRITE (NOUT, 10070, IOSTAT=IDUM)
     &              (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
                endif
                write (nout, 10080)
              ELSE
                if (mxnam .gt. 24) then
                  WRITE (*, 10090, IOSTAT=IDUM)
     &              (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
                else 
                  WRITE (*, 10070, IOSTAT=IDUM)
     &              (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
                endif
                write (*, 10080)
              END IF
            end if
          END IF
          
        IF (DOFAC) THEN
           IF (NDNPS(INPS) .GT. 0) THEN
              IS = IXDNPS(INPS)
              IE = IS + NDNPS(INPS) - 1
C ... See if all values are the same
              val = facnps(is)
              allsam = .TRUE.
              do 50 i=is+1, ie
                 if (facnps(i) .ne. val) then
                    allsam = .FALSE.
                    go to 90
                 end if
 50           continue
 90           continue

              if (allsam) then
                 IF (NOUT .GT. 0) THEN
                    WRITE (NOUT, 10055, IOSTAT=IDUM) VAL
                 ELSE
                    WRITE (*, 10055, IOSTAT=IDUM) VAL
                 END IF
              else
                 IF (NOUT .GT. 0) THEN
                    WRITE (NOUT, 10050, IOSTAT=IDUM)
     &                   (FACNPS(I), I=IS,IE)
                 ELSE
                    WRITE (*, 10050, IOSTAT=IDUM)
     &                   (FACNPS(I), I=IS,IE)
                 END IF
              end if
           ELSE
              IF (NOUT .GT. 0) THEN
                 WRITE (NOUT, 10060, IOSTAT=IDUM)
              ELSE
                 WRITE (*, 10060, IOSTAT=IDUM)
              END IF
           END IF
        END IF
 100  CONTINUE
      
      RETURN

10005 FORMAT (1X, 'Nodal ids are Global')
10020 FORMAT (/, 1X, 'NODAL POINT SETS', :, ' - ', A)
10030 FORMAT (1X, 'Set', I10, 1X, A, ':',
     &  I10, ' nodes', 1X, A, '     name = "',A,'"')
10040 FORMAT ((1X, 8I10))
10050 FORMAT ((1X, 6 (1X, 1pE11.4)))
10055 FORMAT (10x, 'All distribution factors are equal to ', 1pe11.4)
10060 FORMAT (10x, 'Distribution factors not stored in file.')
10070 FORMAT ((2x,4(2X, A)))
10090 FORMAT ((2x,3(2X, A)))
10080 FORMAT (1X)
      END
      
