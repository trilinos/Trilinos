C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE PRELEM (OPTION, NOUT,
     &   NELBLK, NUMEL, NLISEL, LISEL, LENE,
     &   NVAREL, LISEV, NAMEEV, ISEVOK, VAREL, MAPEL)
C=======================================================================

C   --*** PRELEM *** (BLOT) Display current database element variables
C   --   Written by Amy Gilkey - revised 04/20/88
C   --
C   --PRELEM displays the element data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements
C   --   NLISEL - IN - the number of selected elements by block
C   --   LISEL - IN - the indices of the selected elements by block
C   --   LENE - IN - the cumulative element counts by element block
C   --   NVAREL - IN - the number of element variables
C   --   LISEV - IN - the indices of the selected element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VAREL - IN - the selected element variables for the time step

      CHARACTER*(*) OPTION
      INTEGER NLISEL(0:*)
      INTEGER LISEL(0:*)
      INTEGER LENE(0:*)
      INTEGER LISEV(0:*)
      CHARACTER*(*) NAMEEV(*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      REAL VAREL(NUMEL,*)
      INTEGER MAPEL(*)

      LOGICAL ISABRT

      CHARACTER*13 CVAL(5)

      if (nout .gt. 0) then
        WRITE (NOUT, 10000)
      else
        WRITE (*, 10000)
      end if
      
      do i=1, lisev(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMEEV(LISEV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMEEV(LISEV(I))
        END IF
      end do

      DO 130 IELB = 1, NELBLK
         IF (NLISEL(IELB) .GT. 0) THEN

            IX0 = LENE(IELB-1)
            DO 120 IX = 1, NLISEL(IELB)
               IF (ISABRT ()) RETURN
               IEL = LISEL(IX0+IX)

               DO 110 IVAR = 1, LISEV(0), 5
                  MINVAL = IVAR
                  MAXVAL = MIN (LISEV(0), IVAR+5-1)
                  NVAL = MAXVAL - MINVAL + 1
                  DO 100 I = MINVAL, MAXVAL
                     IF (ISEVOK (IELB, LISEV(I))) THEN
                        WRITE (CVAL(I-MINVAL+1), 10020, IOSTAT=IDUM)
     &                     VAREL(IEL,I)
                     ELSE
                        CVAL(I-MINVAL+1) = '-----------'
                     END IF
  100             CONTINUE

                  IF (IVAR .EQ. 1) THEN
                     IF (NOUT .GT. 0) THEN
                       WRITE (NOUT, 10030, IOSTAT=IDUM)
     &                   MAPEL(IEL), (CVAL(I), I=1,NVAL)
                     ELSE
                       WRITE (*, 10030, IOSTAT=IDUM)
     &                   MAPEL(IEL), (CVAL(I), I=1,NVAL)
                     END IF
                  ELSE
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, 10040, IOSTAT=IDUM)
     &                     (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, 10040, IOSTAT=IDUM)
     &                     (CVAL(I), I=1,NVAL)
                     END IF
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
  130 CONTINUE

      RETURN

10000  FORMAT (/, ' Element Time Step Variables (Global Element Ids)')
10010  FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
10020  FORMAT (1PE13.6)
10030  FORMAT (1X, 'Elem', I9, 5 (1X, A13))
10040  FORMAT (1X, 13X, 5 (1X, A13))
      END
