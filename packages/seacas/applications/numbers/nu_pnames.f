C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

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

      NROW = MAX(NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL)

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
      END
