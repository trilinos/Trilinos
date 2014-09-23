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

C $Id: pnames.f,v 1.1 1991/02/21 15:44:51 gdsjaar Exp $
C $Log: pnames.f,v $
C Revision 1.1  1991/02/21 15:44:51  gdsjaar
C Initial revision
C
      SUBROUTINE PNAMES(NAMECO, NAMEBL, NAMEHV, NAMEGV, NAMENV, NAMEEV,
     &               NDIM, NELBLK, NVARHI, NVARGL, NVARNP,NVAREL,COPY)
      CHARACTER*8 NLIST(6), BLANK, NAMECO(*), NAMEBL(*), NAMEHV(*),
     &            NAMEGV(*), NAMENV(*), NAMEEV(*)
      include 'nu_io.blk'
      LOGICAL COPY
C
      DATA BLANK/'        '/
C
C************************************************************************
C
C       G. D. Sjaardema, 1521,  01/30/88
C
C DESCRIPTION: Read and transfer the names found on the data base and 
C       print a formatted echo to SYS$OUTPUT
C
C DUMMY VARIABLES: 
C       NAMECO    CHARACTER     Names of coordinates
C       NAMEBL    CHARACTER     Names of element blocks
C       NAMEHV    CHARACTER     Names of history variables
C       NAMEGV    CHARACTER     Names of global variables
C       NAMENV    CHARACTER     Names of nodal variables
C       NAMEEV    CHARACTER     Names of element variables
C       NDIM      INTEGER       Number of spatial dimensions
C       NELBLK    INTEGER       Number of element blocks
C       NVARHI    INTEGER       Number of history variables
C       NVARGL    INTEGER       Number of global variables
C       NVARNP    INTEGER       Number of nodal variables
C       NVAREL    INTEGER       Number of element variables
C       COPY      LOGICAL       TRUE if echo to output data base
C
C COMMON VARIABLES: --NONE--
C
C FILES:
C       UNIT NDB - INPUT, SEQUENTIAL, UNFORMATTED, READONLY
C       UNIT 11 - OUTPUT, SEQUENTIAL, UNFORMATTED
C               - Output database, written iff COPY = .TRUE.
C
C INTRINSICS CALLED:
C       MAX -- Get maximum value of items in list
C
C ROUTINES CALLED: --NONE--
C
C************************************************************************
C
         READ  (NDB,END=2000,ERR=2100) (NAMEHV(I),I=1,NVARHI),
     $              (NAMEGV(I),I=1,NVARGL),
     $              (NAMENV(I),I=1,NVARNP),
     $              (NAMEEV(I),I=1,NVAREL)
C      PRINT 1000
 1000 FORMAT (/T6,'Coordinate',T18,'History',T30,'Global',T42,
     *            'Nodal',T54,'Element',T66,'Block',
     *        /T6,'----------',T18,'-------',T30,'------',T42,
     *            '-----',T54,'-------',T66,'-----')
C
      NROW = MAX(NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL)
C
      IF (.FALSE.) THEN
      DO 10 I=1, NROW+1
         DO 5 J=1,6
            NLIST(J) = BLANK
   5     CONTINUE
C -COORDINATE NAMES
         IF (I .LE. NDIM)   NLIST(1) = NAMECO(I)
C -HISTORY NAMES
         IF (I .LE. NVARHI) NLIST(2) = NAMEHV(I)
C -GLOBAL NAMES
         IF (I .LE. NVARGL) NLIST(3) = NAMEGV(I)
C -NODAL VARIABLE NAMES
         IF (I .LE. NVARNP) NLIST(4) = NAMENV(I)
C -ELEMENT NAMES
         IF (I .LE. NVAREL) NLIST(5) = NAMEEV(I)
C -ELEMENT BLOCK NAMES
         IF (I .LE. NELBLK) NLIST(6) = NAMEBL(I)
C
         PRINT 1500, (NLIST(J),J=1,6)
 1500    FORMAT (T6,A8,T18,A8,T30,A8,T42,A8,T54,A8,T66,A8)
   10 CONTINUE
      END IF
C
       RETURN
C
C END OF FILE OR READ/WRITE ERROR DURING TRANSFER
C
 2000 CONTINUE
      PRINT *, 'End of file during names transfer'
      STOP 'PNAMES'
 2100 CONTINUE
      PRINT *, 'Read error during names transfer'
      STOP 'PNAMES'
 2200 CONTINUE 
      PRINT *, 'Write error during names transfer'
      STOP 'PNAMES'
      END      
