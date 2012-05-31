C $Id: putcrs.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: putcrs.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:13:57  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:13:56  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]PUTCRS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PUTCRS (X, Y, OLDCUR)
C***********************************************************************
C
C     SUBROUTINE PUTCRS = PLACES THE CROSSHAIRS AT THE CURRENT LOCATION
C
C***********************************************************************
C
C
      DIMENSION IDUM(2)
C
      LOGICAL OLDCUR
C
      CHARACTER DUMMY*16
C
C  SELECT DECIMAL MODE
C
      DUMMY = CHAR(27)
      DUMMY(2:4) = 'OR1'
      WRITE(*,*)DUMMY
C
C  PLACE THE CROSSHAIRS AT THE RIGHT LOCATION
C
      CALL MP2PT(1, X, Y, X1, Y1, IDUM)
      IX = INT(X1*4151.)
      IY = INT(Y1*4151.)
      DUMMY(1:1) = CHAR(27)
      DUMMY(2:2) = 'P'
      WRITE(DUMMY(3:8), '(I6)')IX
      DUMMY(9:9) = ','
      WRITE(DUMMY(10:15), '(I6)')IY
      DUMMY(16:16) = ','
      WRITE(*,*)DUMMY
C
C  UNSELECT DECIMAL MODE
C
      DUMMY = CHAR(27)
      DUMMY(2:4) = 'OR0'
      WRITE(*,*)DUMMY
C
      IF(.NOT.OLDCUR)THEN
C
C  ACTIVATE THE CROSSHAIRS
C
         DUMMY = CHAR(27)
         DUMMY(2:3) = 'G1'
         WRITE(*,*)DUMMY
         OLDCUR = .TRUE.
      ENDIF
C
      WRITE(*, '(A)')' '//CHAR(27)//'[2J'
      RETURN
C
      END
