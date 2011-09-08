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
      SUBROUTINE MYRSRV (MYCV, NAME1, NEWLEN, NEWLOC, MYLOC, MYCLOC,
     *   UCLOC, OFFSET, COFFST, VOID, LVOID,
     *   NVOIDS, DICT, DPOINT, LDICT, NNAMES, CHRCOL, CHRNUM,
     *   DEFER, CFILL, CFDATA, MAXSIZ,
     *   LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This routine finds space to service a non-negative space request.
C     If zero space is requested, a valid pointer of 1 will be
C     generated.
C
C***********************************************************************
C
C     MYCV     Internal reference array.
               CHARACTER MYCV(*)
C     NAME1    Name to be inserted in the dictionary
               CHARACTER*8 NAME1
C     NEWLEN   Length of requested storage (character units)
C     NEWLOC   Pointer of new storage (returned)
C     MYLOC    Reference address of internal numeric array
C     MYCLOC   Address of internal character array.
C     UCLOC    Address of user's character array.
C     OFFSET   Offset between internal numeric array and user's
C              numeric array.
C     COFFST   Offset between internal numeric array and user's
C              character array.
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2)
               DIMENSION NVOIDS(2)
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary tables
C     NNAMES   Number of names
               CHARACTER DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3)
               DIMENSION NNAMES(2)
C     CHRCOL   Number of column for character names.
C     CHRNUM   Number of characters per numeric storage unit.
C     DEFER    Flag for deferred mode.
               LOGICAL DEFER
C     CFILL    Flag for character data fill.
C     CFDATA   Data for fill.
               LOGICAL CFILL
               CHARACTER*1 CFDATA
C     MAXSIZ   Dimension of static character array.
C     LASTER   Error return
C
C***********************************************************************
C
      LASTER = SUCESS
      INTLEN = (NEWLEN + CHRNUM - 1) / CHRNUM
C
      IF (NEWLEN .EQ. 0) THEN
C
C        Zero length entry.
C
         NEWLOC = 1 - COFFST / CHRNUM
      ELSE
C
         CALL MXLOOK (INTLEN, VOID(1,CHRCOL,1), CHRCOL*LVOID,
     *      NVOIDS(CHRCOL), VROW, LASTER)
C
         IF (LASTER .EQ. SUCESS) THEN
            NEWLOC = VOID(VROW,1,1)
         ELSE IF (DEFER .AND. CHRCOL .EQ. 1) THEN
C
C           A good void was not found - defer the space request.
C
            NEWLOC = IXLNUM(NEWLOC)
            INTLEN = - INTLEN
            LASTER = SUCESS
C
         ELSE IF (CHRCOL .EQ. 1) THEN
C
C           Get space.
C
            CALL MXGET (MYLOC, INTLEN, VOID, LVOID,
     *         NVOIDS, CHRCOL, LASTER, VROW)
            IF (LASTER .NE. SUCESS) RETURN
            NEWLOC = VOID(VROW,1,1)
C
         ELSE
C
C           CHRCOL .EQ. 2
C
            CALL MYGET (MYCLOC, NEWLEN, VOID, LVOID,
     *         NVOIDS, CHRCOL, MAXSIZ, LASTER, VROW)
            IF (LASTER .NE. SUCESS) RETURN
            NEWLOC = VOID(VROW,2,1)
C
         END IF
      END IF
C
C     Update dictionary.
C
      CALL MYNSRT (NAME1, NEWLOC, INTLEN, NEWLEN, DICT, DPOINT, LDICT,
     *   NNAMES, CHRCOL, LASTER)
      IF (LASTER .EQ. WRTYPE) LASTER = BDNAME
      IF (LASTER .NE. SUCESS) RETURN
C
      IF (INTLEN .GT. 0) THEN
C
C        Data fill pattern.
C
         IF (CFILL) THEN
            TLOC = (VOID(VROW,CHRCOL,1) - 1) * CHRNUM + 1 + COFFST
     *         + UCLOC - MYCLOC
            DO 100 I = TLOC, TLOC + NEWLEN - 1
               MYCV(I) = CFDATA
  100       CONTINUE
         END IF
C
C        Update void table.
C
         VOID(VROW,CHRCOL,1) = VOID(VROW,CHRCOL,1) + INTLEN
         VOID(VROW,CHRCOL,2) = VOID(VROW,CHRCOL,2) - INTLEN
         CALL VTABLE (1, 0, VOID(1,CHRCOL,1), LVOID, NVOIDS(CHRCOL),
     *      CHRCOL, LASTER)
         NEWLOC = (NEWLOC - 1) * CHRNUM + 1 + COFFST
      ELSE
         NEWLOC = - UCLOC
      END IF
C
      RETURN
      END
