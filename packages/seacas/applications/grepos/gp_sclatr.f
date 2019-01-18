C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE SCLATR (NELBLK, IDELB, NUMATR, IDATR, DOALLA,
     *  IDBLK, DOALLB, ASCALE, ATRSCL, VERBOS, SCLSET)
C=======================================================================

      INTEGER IDELB(*)
      INTEGER NUMATR(*)
      LOGICAL DOALLA
      LOGICAL DOALLB
      REAL    ATRSCL(2,*)
C ... Row 1 == value to set attribute to, Row 2 == value to scale by.
      LOGICAL VERBOS
      LOGICAL SCLSET
C ... SCLSET == .TRUE. then do scale, .FALSE. then set to value      

      IAT = 0
      DO 110 IELB = 1, NELBLK
        IF (IDELB(IELB) .EQ. IDBLK .OR. DOALLB) THEN
          DO 100 IATR = 1, NUMATR(IELB)
            IF (IDATR .EQ. IATR .OR. DOALLA) THEN
              IF (SCLSET) THEN
                ATRSCL(1,IAT+IATR) = 0.0
                ATRSCL(2,IAT+IATR) = ASCALE
                if (verbos) then
                  WRITE (*, 10000) IATR, IDELB(IELB), ' Scaled by ',
     *              ASCALE
                end if
              ELSE
                ATRSCL(1,IAT+IATR) = ASCALE
                ATRSCL(2,IAT+IATR) = 1.0
                if (verbos) then
                  WRITE (*, 10000) IATR, IDELB(IELB), ' Set to ', ASCALE
                end if
              END IF
              
            END IF
 100      CONTINUE
        END IF
        IAT = IAT + NUMATR(IELB)
 110  CONTINUE
      RETURN

10000 FORMAT (1X, 'Attribute', I3, ' of Block', I5,A, 1PE10.3)
      END

      SUBROUTINE CNTATR (NELBLK, NUMATR, IAT)
      INTEGER NUMATR(*)
      IAT = 0
      DO 110 IELB = 1, NELBLK
        IAT = IAT + NUMATR(IELB)
 110  CONTINUE
      RETURN
      END
