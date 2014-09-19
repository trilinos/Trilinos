C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: selbox.f,v 1.2 1998/03/22 05:34:43 gdsjaar Exp $
C $Log: selbox.f,v $
C Revision 1.2  1998/03/22 05:34:43  gdsjaar
C General cleanp of unused variables. Reordered DATA statements in
C command.f so would compile with f2c.
C
C Revision 1.1.1.1  1991/02/21 15:45:24  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:23  gdsjaar
c Initial revision
c
      SUBROUTINE SELBOX (COORD, NUMNP, NDIM, P1, SELECT, NODEL)
      DIMENSION COORD (NUMNP,*), P1(*)
      LOGICAL SELECT(*)
      CHARACTER*8 NODEL
      CHARACTER*80 STRTMP
      INTEGER LENSTR
C
      CALL INILOG (NUMNP, .FALSE., SELECT)
      INUM = 0
      IF (NDIM .EQ. 2) THEN
      DO 10 I=1, NUMNP
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            IF (X0 .GE. P1(1) .AND. X0 .LE. P1(2) .AND.
     *          Y0 .GE. P1(3) .AND. Y0 .LE. P1(4)) THEN
               SELECT(I) = .TRUE.
               INUM = INUM + 1
         END IF
   10 CONTINUE
      ELSE IF (NDIM .EQ. 3) THEN
      DO 20 I=1, NUMNP
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            Z0 = COORD(I,3)
            IF (X0 .GE. P1(1) .AND. X0 .LE. P1(2) .AND.
     *          Y0 .GE. P1(3) .AND. Y0 .LE. P1(4) .AND. 
     *          Z0 .GE. P1(5) .AND. Z0 .LE. P1(6)) THEN
               SELECT(I) = .TRUE.
               INUM = INUM + 1
         END IF
   20 CONTINUE
      ELSE
            CALL PRTERR ('PROGRAM', 'Illegal dimension in SELBOX')
      END IF

      IF (INUM .EQ. 0) THEN
         CALL PRTERR ('WARNING', 
     *        'No '// NODEL(:LENSTR(NODEL))//' found in range.')
      ELSE
         WRITE (STRTMP, 100) INUM, NODEL
  100    FORMAT (I10,' ',A,' Selected.')
         CALL SQZSTR (STRTMP, LTMP)
         CALL PRTERR ('CMDSPEC', STRTMP(:LTMP))
      END IF
      RETURN
      END
