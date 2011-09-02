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
      SUBROUTINE PRNSV (NOUT, NSTEP, NUMNPS, LISNPS, LNPSNL,
     &     IDNPS, NNNPS, IXNNPS, LTNNPS, NAME, 
     $     NVAR,  LISVAR, NAMEV, ISVOK, VARS, nvardm,
     $     MAPNO, DOMAP)
C=======================================================================

C     --*** PRNPS *** (GROPE) Display database nodal point set
C     --
C     --PRNPS displays the nodal point sets.
C     --
C     --Parameters:
C     --   NOUT - IN - the output file, <=0 for standard
C     --   NUMNPS - IN - the number of nodal point sets
C     --   LISNPS - IN - the indices of the selected nodal point sets
C     --   LNPSNL - IN - the number of nodes for all sets
C     --   IDNPS - IN - the nodal point set ID for each set
C     --   NNNPS - IN - the number of nodes for each set
C     --   NDNPS - IN - the number of distribution factors for each set
C     --   IXNNPS - IN - the index of the first node for each set
C     --   IXDNPS - IN - the index of the first dist factor for each set
C     --   LTNNPS - IN - the nodes for all sets
C     --   FACNPS - IN - the distribution factors for all sets

      include 'params.blk'
      include 'dbase.blk'
      INTEGER LISNPS(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      CHARACTER*(*) NAME(*)

      INTEGER LISVAR(0:*)
      CHARACTER*(*) NAMEV(*)
      INTEGER ISVOK(nvardm,*)
      REAL VARS(LNPSNL, *)

      INTEGER MAPNO(*)
      LOGICAL DOMAP

      CHARACTER*32 CVAL(6)
      CHARACTER*40 FMT20,FMT30, FMT40
      INTEGER PRTLEN
      INTEGER GETPRC

      CHARACTER*20 STRA, STRB

C ... See if need to read the data
      if (nstep .ne. nstepns) then
         nstepns = nstep
         DO 20 IX = 1, LISNPS(0)
            INPS = LISNPS(IX)
            ID   = idnps(inps)
            IS = IXNNPS(INPS)
            DO 10 IVAR = 1, LISVAR(0)
               IF (ISVOK (LISVAR(IVAR),INPS) .NE. 0) THEN
                  call exgnsv(ndb, nstep, ivar, id, nnnps(inps),
     $                 vars(is,ivar), ierr)
               end if
 10         continue
            ioff = ioff + nnnps(inps)
 20      continue
      end if

      PRTLEN = GETPRC() + 7
      WRITE(FMT20,2000) PRTLEN, PRTLEN-7
      WRITE(FMT30,3000) PRTLEN
      WRITE(FMT40,4000) PRTLEN

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      if (domap) then
         if (nout .gt. 0) then
            write (nout, 10005)
         else
            write (*, 10005)
         end if
      end if

      do 90 i=1, lisvar(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMEV(LISVAR(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMEV(LISVAR(I))
        END IF
 90   continue


      WRITE (STRA, 10001, IOSTAT=IDUM) NUMNPS
10001 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10002, IOSTAT=IDUM) LNPSNL
10002 FORMAT ('(index=', I10, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      DO 200 IX = 1, LISNPS(0)
         INPS = LISNPS(IX)
         WRITE (STRA, 10001, IOSTAT=IDUM) INPS
         CALL PCKSTR (1, STRA)
         WRITE (STRB, 10002, IOSTAT=IDUM) IXNNPS(INPS)
         CALL PCKSTR (1, STRB)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &           IDNPS(INPS), STRA(:LSTRA),
     &           NNNPS(INPS), STRB(:LSTRB),
     $           NAME(INPS)(:LENSTR(NAME(INPS)))
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &           IDNPS(INPS), STRA(:LSTRA),
     &           NNNPS(INPS), STRB(:LSTRB),
     $           NAME(INPS)(:LENSTR(NAME(INPS)))
         END IF

         IF (NNNPS(INPS) .GT. 0) THEN
            IS = IXNNPS(INPS)
            IE = IS + NNNPS(INPS) - 1
            do 130 IN = IS,IE
               if (domap) then
                  ID = MAPNO(LTNNPS(IN))
               else
                  ID = LTNNPS(IN)
               end if
               
               DO 110 IVAR = 1, LISVAR(0), 5
                  MINVAL = IVAR
                  MAXVAL = MIN (LISVAR(0), IVAR+5-1)
                  NVAL = MAXVAL - MINVAL + 1
                  DO 100 I = MINVAL, MAXVAL
                     IF (ISVOK (LISVAR(I),INPS) .NE. 0) THEN
                        WRITE (CVAL(I-MINVAL+1), FMT20, IOSTAT=IDUM)
     &                       VARS(IN, LISVAR(I))
                     ELSE
                        CVAL(I-MINVAL+1) = '-----------'
                     END IF
 100              CONTINUE

                  IF (IVAR .LE. 1) THEN
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT30, IOSTAT=IDUM)
     &                       ID, (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT30, IOSTAT=IDUM)
     &                       ID, (CVAL(I), I=1,NVAL)
                     END IF
                  ELSE
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT40, IOSTAT=IDUM)
     &                       (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT40, IOSTAT=IDUM)
     &                       (CVAL(I), I=1,NVAL)
                     END IF
                  END IF
 110           CONTINUE
 130        CONTINUE
         end if
 200  CONTINUE
      
      RETURN

10000 FORMAT (/, 1X, 'NODESET TIME STEP VARIABLES')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
 2000 FORMAT('(1PE',I2.2,'.',I2.2,')')
 3000 FORMAT('(1X, ''Node'', I10, 5(2X,A',I2,'))')
 4000 FORMAT('(15X, 5 (2X, A',I2,'))')

10005 FORMAT (1X, 'Nodal ids are Global')
10020 FORMAT (/, 1X, 'NODAL POINT SETS', :, ' - ', A)
10030 FORMAT (/,1X, 'Set', I10, 1X, A, ':',
     &     I10, ' nodes', 1X, A, ' name = "',A,'"')
10050 FORMAT ((1X, 6 (1X, 1pE11.4)))
      END
      
