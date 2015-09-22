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
      SUBROUTINE MXRSRV (MYV, NAME1, NEWLEN, NEWLOC, MYLOC, OFFSET,
     *   VOID, LVOID,
     *   NVOIDS, DICT, DPOINT, LDICT, NNAMES, CHRCOL,
     *   DEFER, FILL, FDATA,
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
C     MYV      Internal reference array.
               DIMENSION MYV(*)
C     NAME1    Name to be inserted in the dictionary
               CHARACTER*8 NAME1
C     NEWLEN   Length of requested storage
C     NEWLOC   Pointer of new storage (returned)
C     MYLOC    Reference address of internal array
C     OFFSET   Offset between internal array and user's array
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary tables
C     NNAMES   Number of names
               CHARACTER DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     DEFER    Flag for deferred mode.
               LOGICAL DEFER
C     FILL     Flag for data fill.
C     FDATA    Data for fill.
               LOGICAL FILL
C     LASTER   Error return
C
C***********************************************************************
C
      LASTER = SUCESS
      MYLEN = NEWLEN
C
      IF (NEWLEN .EQ. 0) THEN
C
C        Zero length entry.
C
         NEWLOC = 1 - OFFSET
C
      ELSE IF (DEFER) THEN
C
         CALL MXLOOK (MYLEN, VOID, CHRCOL*LVOID, NVOIDS(1),
     *      VROW, LASTER)
C
         IF (LASTER .EQ. SUCESS) THEN
            NEWLOC = VOID(VROW,1,1)
         ELSE IF (LASTER .EQ. NOGET) THEN
C
C           A good void was not found - defer the space request.
C
            NEWLOC = IXLNUM(NEWLOC)
            MYLEN = - NEWLEN
            LASTER = SUCESS
C
         END IF
C
      ELSE
C
C        Get space.
C
         CALL MXGET (MYLOC, MYLEN, VOID, LVOID, NVOIDS,
     *      CHRCOL, LASTER, VROW)
         IF (LASTER .NE. SUCESS) RETURN
C
         NEWLOC = VOID(VROW,1,1)
C
      END IF
C
C     Update dictionary.
C
      CALL MXNSRT (NAME1, NEWLOC, MYLEN, DICT, DPOINT, LDICT,
     *   NNAMES, CHRCOL, LASTER)
      IF (LASTER .EQ. WRTYPE) LASTER = BDNAME
      IF (LASTER .NE. SUCESS) RETURN
C
      IF (MYLEN .GT. 0) THEN
C
C        Data fill pattern.
C
         IF (FILL) THEN
            DO 100 I = VOID(VROW,1,1), VOID(VROW,1,1)+MYLEN-1-7,8
               MYV(I+0) = FDATA
               MYV(I+1) = FDATA
               MYV(I+2) = FDATA
               MYV(I+3) = FDATA
               MYV(I+4) = FDATA
               MYV(I+5) = FDATA
               MYV(I+6) = FDATA
               MYV(I+7) = FDATA
  100       CONTINUE
            do 110 J = I, VOID(VROW,1,1)+MYLEN-1
              MYV(J) = FDATA
 110        continue
         END IF
C
C        Update void table.
C
         VOID(VROW,1,1) = VOID(VROW,1,1) + MYLEN
         VOID(VROW,1,2) = VOID(VROW,1,2) - MYLEN
         CALL VTABLE (1, 0, VOID, LVOID, NVOIDS(1), CHRCOL, LASTER)
         NEWLOC = NEWLOC + OFFSET
C
      ELSE IF (MYLEN .LT. 0) THEN
C
         NEWLOC = OFFSET - MYLOC
C
      END IF
C
      RETURN
      END
