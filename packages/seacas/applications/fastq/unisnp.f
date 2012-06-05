C $Id: unisnp.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: unisnp.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:17:32  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:17:30  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]UNISNP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, XMIN, XMAX, STEP)
C***********************************************************************
C
C  SUBROUTINE UNISNP = GENERATES UNIFORM SNAP GRID
C
C***********************************************************************
C
C  VARIABLES USED:
C     MSNAP  = DIMENSION OV SNAP ARRAYS
C     SNAPDX = THE SNAP GRID VALUES ARRAY  (X AND Y)
C     NSNAP  = THE NUMBER OF SNAP GRID VALUES IN X AND Y
C     INDEX  = 1 FOR X VALUES,  2 FOR Y VALUES
C
C***********************************************************************
C
      DIMENSION SNAPDX (2, MSNAP), NSNAP (2)
C
      LOGICAL ERR
C
      CHARACTER*1 AXIS (2)
C
      DATA AXIS /'X', 'Y'/
C
C  DEFINE THE GRID
C
      IF (STEP.EQ.0.)RETURN
      ILOOP =  INT(((XMAX - XMIN) / STEP) + 2)
      XGRID = XMIN
      DO 100 I = 1, ILOOP
         CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, XGRID, ERR)
         IF (ERR)THEN
            WRITE (*, 10000)AXIS (INDEX),  XGRID - STEP
            RETURN
         ENDIF
         XGRID = XGRID + STEP
         IF (XGRID.GE. (STEP + XMAX))RETURN
  100 CONTINUE
C
      RETURN
C
10000 FORMAT (' THE LAST SUCESSFULL ', A1, ' GRID INPUT WAS: ', G14.7)
C
      END
