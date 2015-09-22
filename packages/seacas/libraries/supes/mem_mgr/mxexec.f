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
      SUBROUTINE MXEXEC (MYV, MYCV, MYLOC, MYCLOC, UCLOC, COFFST,
     *   OFFSET, DPOINT, LDICT, NNAMES,
     *   VOID, LVOID, NVOIDS, FILL, FDATA, CFILL, CFDATA, CHRNUM,
     *   CHRCOL, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This routine satisfies deferred memory requests.
C     It will service both numeric and character
C     deferred requests in the mixed mode.  In the nonmixed mode,
C     character memory is not deferred.
C
C***********************************************************************
C
C     MYV      Internal reference vector.
               DIMENSION MYV(*)
C     MYCV     Internal reference array.
               CHARACTER MYCV(*)
C     MYLOC    Address of internal reference vector.
C     MYCLOC   Address of internal character array.
C     UCLOC    Address of user's character array.
C     COFFST   Offset between internal numeric array and user's
C              character array.
C     OFFSET   Offset address between internal and user arrays.
C     DPOINT   Dictionary pointer table.
C     LDICT    Dimension of dictionary table.
C     NNAMES   Number of names in dictionary.
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     VOID     Void table.
C     LVOID    Dimension of void table.
C     NVOIDS   Number of voids.
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     FILL     Flag for data fill.
C     FDATA    Data for fill.
               LOGICAL FILL
C     CFILL    Flag for character data fill.
C     CFDATA   Data for fill.
               LOGICAL CFILL
               CHARACTER*1 CFDATA
C     CHRNUM   Number of characters per numeric storage unit
C     CHRCOL   Column for character tables.
C     LASTER   Error return.
C
C***********************************************************************
C
      LASTER = SUCESS
      MEM = 0
      DO 100 IDICT = 1, NNAMES(1)
         MEM = MEM + MIN(0, DPOINT(IDICT,1,2))
  100 CONTINUE
      IF (MEM .EQ. 0) RETURN
      MEM = - MEM
C
C
      CALL MXGET (MYLOC, MEM, VOID, LVOID, NVOIDS,
     *   CHRCOL, LASTER, VROW)
      IF (LASTER .NE. SUCESS) RETURN
C
C     Now satisfy all the deferred requests.
C
      DO 130 IDICT = 1, NNAMES(1)
         IF (DPOINT(IDICT,1,2) .LT. 0) THEN
            IF (DPOINT(IDICT,1,3) .EQ. -1) THEN
               MYV(DPOINT(IDICT,1,1) - MYLOC + 1) = 
     *            VOID(VROW,1,1) + OFFSET
            ELSE
               MYV(DPOINT(IDICT,1,1) - MYLOC + 1) = 
     *            (VOID(VROW,1,1) - 1) * CHRNUM + 1 + COFFST
            END IF
            VOID(VROW,1,2) = VOID(VROW,1,2) + DPOINT(IDICT,1,2)
            DPOINT(IDICT,1,2) = - DPOINT(IDICT,1,2)
            DPOINT(IDICT,1,1) = VOID(VROW,1,1)
            VOID(VROW,1,1) = VOID(VROW,1,1) + DPOINT(IDICT,1,2)
C
C           Perform data fill if appropriate.
C
            IF (FILL .AND. DPOINT(IDICT,1,3) .EQ. -1) THEN
               DO 110 I = DPOINT(IDICT,1,1),
     *            DPOINT(IDICT,1,1)+DPOINT(IDICT,1,2)-1
                  MYV(I) = FDATA
  110          CONTINUE
            ELSE IF (CFILL .AND. DPOINT(IDICT,1,3) .GT. 0) THEN
               TLOC = (DPOINT(IDICT,1,1) - 1) * CHRNUM + 1 + COFFST
     *            + UCLOC - MYCLOC
               DO 120 I = TLOC, TLOC + DPOINT(IDICT,1,3) - 1
                  MYCV(I) = CFDATA
  120          CONTINUE
            END IF
C
         END IF
  130 CONTINUE
C
      CALL VTABLE (0, 0, VOID, LVOID, NVOIDS(1), CHRCOL, LASTER)
      RETURN
      END
