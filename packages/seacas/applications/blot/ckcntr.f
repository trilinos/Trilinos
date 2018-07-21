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
      SUBROUTINE CKCNTR (OK)
C=======================================================================

C   --*** CKCNTR *** (DETOUR) Check the contour values
C   --   Written by Amy Gilkey - revised 07/10/87
C   --
C   --CKCNTR checks that all specified contour values increase or decrease.
C   --
C   --Parameters:
C   --   OK - OUT - true iff the contour values are consistent
C   --
C   --Common Variables:
C   --   Uses CINTOK, LINCON, NCNTR, CINTV of /CNTR/

      include 'cntr.blk'
C     FLAG FOR EXACT CONTOUR VALUES FOR EACH PLOT
C     COMMON /CNTR/   CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC,
C    &   CINTV(256), NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX
C     LOGICAL CINTOK, LINCON, NOCMIN, NOCMAX


      LOGICAL OK

      CHARACTER*80 ERRSTR

      OK = .TRUE.

      IF (CINTOK) THEN
         IF (LINCON) THEN
            NC = NCNTR
         ELSE
            NC = NCNTR+1
         END IF

         IF (CINTV(1) .LE. CINTV(2)) THEN
            DO 100 I = 2, NC
               IF (CINTV(I-1) .GE. CINTV(I)) THEN
                  WRITE (ERRSTR, 10000)
     &               'Contour interval ', I-1, ' >= interval ', I
10000             FORMAT (A, I5, A, I5)
                  CALL SQZSTR (ERRSTR, LSTR)
                  CALL PRTERR ('CMDWARN', ERRSTR(:LSTR))
                  OK = .FALSE.
               END IF
  100       CONTINUE
         ELSE
            DO 110 I = 2, NC
               IF (CINTV(I-1) .LE. CINTV(I)) THEN
                  WRITE (ERRSTR, 10000)
     &               'Contour interval ', I-1, ' <= interval ', I
                  CALL SQZSTR (ERRSTR, LSTR)
                  CALL PRTERR ('CMDWARN', ERRSTR(:LSTR))
                  OK = .FALSE.
               END IF
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
