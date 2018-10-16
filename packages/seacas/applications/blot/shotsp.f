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

C $Log: shotsp.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:12:21  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:57:24  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SHOTSP (WHONLY, TMIN, TMAX, DELT, NINTV, NPTIMS)
C=======================================================================

C   --*** SHOTSP *** (TIMSEL) Display time step parameters
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --SHOTSP displays the time step selection parameters.
C   --
C   --Parameters:
C   --   WHONLY - IN - true iff only whole times may be selected
C   --   TMIN - IN - the minimum selected time
C   --   TMAX - IN - the maximum selected time
C   --   DELT - IN - the selected times interval (<0 = selected times)
C   --   NINTV - IN - negative for zero interval
C   --   NPTIMS - IN - the number of selected times

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      LOGICAL WHONLY
      REAL TMIN, TMAX, DELT
      INTEGER NINTV
      INTEGER NPTIMS

      CHARACTER*80 STRING
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      IF (DELT .GT. 0) THEN
         IF (NINTV .NE. 0) THEN
            RNUM(1) = TMIN
            RNUM(2) = TMAX
            CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
            IF (NINTV .GT. 0) THEN
               IF (WHONLY) THEN
                  WRITE (STRING, 10000) ' whole ',
     &               (RSTR(I)(:LSTR), I=1,2), NINTV, 'delta'
               ELSE
                  WRITE (STRING, 10000) ' ',
     &               (RSTR(I)(:LSTR), I=1,2), NINTV, 'delta'
               END IF
            ELSE
               IF (WHONLY) THEN
                  WRITE (STRING, 10000) ' whole ',
     &               (RSTR(I)(:LSTR), I=1,2), -NINTV, 'zero'
               ELSE
                  WRITE (STRING, 10000) ' ',
     &               (RSTR(I)(:LSTR), I=1,2), -NINTV, 'zero'
               END IF
            END IF
10000        FORMAT ('Select', A, 'times ', A, ' to ', A, ' in ', I6,
     &         ' intervals with ', A, ' offset')
         ELSE
            RNUM(1) = TMIN
            RNUM(2) = TMAX
            RNUM(3) = DELT
            CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
            IF (WHONLY) THEN
               WRITE (STRING, 10010) ' whole ', (RSTR(I)(:LSTR), I=1,3)
            ELSE
               WRITE (STRING, 10010) ' ', (RSTR(I)(:LSTR), I=1,3)
            END IF
10010        FORMAT ('Select', A, 'times ', A, ' to ', A, ' by ', A)
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10030) STRING(:LSTR)

      ELSE IF (DELT .EQ. 0) THEN
         RNUM(1) = TMIN
         RNUM(2) = TMAX
         CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
         IF (WHONLY) THEN
            WRITE (STRING, 10020) ' whole ', (RSTR(I)(:LSTR), I=1,2)
         ELSE
            WRITE (STRING, 10020) ' ', (RSTR(I)(:LSTR), I=1,2)
         END IF
10020     FORMAT ('Select all', A, 'times from ', A, ' to ', A)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10030) STRING(:LSTR)

      ELSE
         IF (WHONLY) THEN
            WRITE (*, 10030) 'Select specified whole times'
         ELSE
            WRITE (*, 10030) 'Select specified times'
         END IF
      END IF

      CALL INTSTR (1, 0, NPTIMS, STRING, LSTR)
      WRITE (*, 10030) '   Number of selected times = ', STRING(:LSTR)

      RETURN
10030  FORMAT (1X, 5A)
      END
