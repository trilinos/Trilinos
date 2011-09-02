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
      SUBROUTINE PRSSV (NOUT, NSTEP, NUMESS, LISESS, LESSEL,
     &     IDESS, NEESS, IXEESS, LTEESS, LTSESS, NAME, 
     $     NVAR,  LISVAR, NAMEV, ISVOK, VARS, nvardm,
     $     MAPEL, DOMAP)
C=======================================================================

C     --*** PRSSV *** (GROPE) Display database nodal point set
C     --
C     --PRSSV displays the sides set vars
C     --
C     --Parameters:
C     --   NOUT - IN - the output file, <=0 for standard
C     --   NUMESS - IN - the number of nodal point sets
C     --   LISESS - IN - the indices of the selected side sets
C     --   LESSEL - IN - the number of elements for all sets
C     --   IDESS - IN - the set ID for each set
C     --   NEESS - IN - the number of elements for each set
C     --   IXEESS - IN - the index of the first element for each set
C     --   LTEESS - IN - the elements for all sets
C     --   LTSESS - IN - the element sides for all sets
      include 'params.blk'
      include 'dbase.blk'
      INTEGER LISESS(0:*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      CHARACTER*(*) NAME(*)

      INTEGER LISVAR(0:*)
      CHARACTER*(*) NAMEV(*)
      INTEGER ISVOK(nvardm,*)
      REAL VARS(LESSEL, *)

      INTEGER MAPEL(*)
      LOGICAL DOMAP

      CHARACTER*32 CVAL(6)
      CHARACTER*40 FMT20,FMT30, FMT40
      INTEGER PRTLEN
      INTEGER GETPRC

      CHARACTER*20 STRA, STRB


C ... See if need to read the data
      if (nstep .ne. nstepss) then
         nstepss = nstep
         DO 20 IX = 1, LISESS(0)
            IESS = LISESS(IX)
            IS = IXEESS(IESS)
            IE = IS + NEESS(IESS) - 1

            ID   = idess(iess)
            DO 10 IVAR = 1, LISVAR(0)
               IF (ISVOK (LISVAR(IVAR),IESS) .NE. 0) THEN
                  call exgssv(ndb, nstep, ivar, id, neess(iess),
     $                 vars(is,ivar), ierr)
               end if
 10         continue
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


      WRITE (STRA, 10001, IOSTAT=IDUM) NUMESS
10001 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10002, IOSTAT=IDUM) LESSEL
10002 FORMAT ('(index=', I10, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)
      
      DO 200 IX = 1, LISESS(0)
         IESS = LISESS(IX)
         WRITE (STRA, 10001, IOSTAT=IDUM) IESS
         CALL PCKSTR (1, STRA)
         WRITE (STRB, 10002, IOSTAT=IDUM) IXEESS(IESS)
         CALL PCKSTR (1, STRB)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &           IDESS(IESS), STRA(:LSTRA),
     &           NEESS(IESS), STRB(:LSTRB),
     $           NAME(IESS)(:LENSTR(NAME(IESS)))
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &           IDESS(IESS), STRA(:LSTRA),
     &           NEESS(IESS), STRB(:LSTRB),
     $           NAME(IESS)(:LENSTR(NAME(IESS)))
         END IF

         IF (NEESS(IESS) .GT. 0) THEN
            IS = IXEESS(IESS)
            IE = IS + NEESS(IESS) - 1
            do 130 IN = IS,IE
               if (domap) then
                  ID = MAPEL(LTEESS(IN))
               else
                  ID = LTEESS(IN)
               end if
               
               DO 110 IVAR = 1, LISVAR(0), 5
                  MINVAL = IVAR
                  MAXVAL = MIN (LISVAR(0), IVAR+5-1)
                  NVAL = MAXVAL - MINVAL + 1
                  DO 100 I = MINVAL, MAXVAL
                     IF (ISVOK (LISVAR(I),IESS) .NE. 0) THEN
                        WRITE (CVAL(I-MINVAL+1), FMT20, IOSTAT=IDUM)
     &                       VARS(IN, LISVAR(I))
                     ELSE
                        CVAL(I-MINVAL+1) = '-----------'
                     END IF
 100              CONTINUE

                  IF (IVAR .LE. 1) THEN
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT30, IOSTAT=IDUM)
     &                       ID, LTSESS(IN), (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT30, IOSTAT=IDUM)
     &                       ID, LTSESS(IN), (CVAL(I), I=1,NVAL)
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

10000 FORMAT (/, 1X, 'SIDESET TIME STEP VARIABLES')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
 2000 FORMAT('(1PE',I2.2,'.',I2.2,')')
 3000 FORMAT('(1X, ''Element'', I10,''.'',I1, 5(2X,A',I2,'))')
 4000 FORMAT('(20X, 5 (2X, A',I2,'))')

10005 FORMAT (1X, 'Element ids are Global')
10020 FORMAT (/, 1X, 'SIDE SETS', :, ' - ', A)
10030 FORMAT (/,1X, 'Set', I10, 1X, A, ':',
     &     I10, ' elements', 1X, A, ' name = "',A,'"')
10050 FORMAT ((1X, 6 (1X, 1pE11.4)))
      END
      
