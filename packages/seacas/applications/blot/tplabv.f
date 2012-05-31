C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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
      SUBROUTINE TPLABV (NPLT, IVAR, NAME, NE, LABSTR, MAPEL, MAPND)
C=======================================================================

C   --*** TPLABV *** (TPLOT) Get a plot label
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --TPLABV returns a plot label.
C   --
C   --Parameters:
C   --   NPLT - IN - the label type:
C   --      -1 = name and number (long form)
C   --       0 = name and number (short form)
C   --   IVAR - IN - the variable index which determines type
C   --   NAME - IN - the variable name
C   --   NE - IN - the variable number
C   --   LABSTR - OUT - the label string
C   --
C   --Common Variables:
C   --   Uses NVARNP, NVAREL of /DBNUMS/

      CHARACTER*(*) NAME
      CHARACTER*(*) LABSTR
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER TYP

      CALL DBVTYP (IVAR, TYP, IDUM)

      IF (NPLT .LE. -1) THEN
         IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME
         ELSE IF (TYP .EQ. 'N') THEN
           WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'at NODE',
     *       MAPND(NE)
         ELSE IF (TYP .EQ. 'E') THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'at ELEMENT',
     *       MAPEL(NE)
         END IF

      ELSE

         IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME
         ELSE IF (TYP .EQ. 'N') THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'NODE',
     *       MAPND(NE)
         ELSE IF (TYP .EQ. 'E') THEN
           WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'ELEM',
     *       MAPEL(NE)
         END IF
      END IF

      CALL SQZSTR (LABSTR, LSTR)

      RETURN
10000  FORMAT (A, :, ' ', A, I6)
      END
