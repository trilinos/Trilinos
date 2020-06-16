C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltarr.f,v 1.1 1993/07/16 16:47:39 gdsjaar Exp $
C $Log: pltarr.f,v $
C Revision 1.1  1993/07/16 16:47:39  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTARR(XTAIL,YTAIL,XHEAD,YHEAD,THETA,ARRLEN)
      REAL XTAIL,YTAIL
      REAL XHEAD,YHEAD
      REAL THETA
      REAL ARRLEN
      DATA SMALL/1.E-3/

      IF (ABS(XTAIL-XHEAD).LE.SMALL .AND.
     *    ABS(YTAIL-YHEAD).LE.SMALL) THEN
         RETURN

      END IF

      CALL PLTVCT(1,XTAIL,YTAIL,XHEAD,YHEAD)
      IF (ARRLEN.LE.0.) THEN
         RETURN

      END IF

      DX = XHEAD - XTAIL
      DY = YHEAD - YTAIL
      ALPHA = ATAN2(DY,DX)
      A = -ARRLEN*COS(ALPHA)*COS(THETA) + ARRLEN*SIN(ALPHA)*SIN(THETA) +
     *    XHEAD
      B = -ARRLEN*COS(ALPHA)*SIN(THETA) - ARRLEN*SIN(ALPHA)*COS(THETA) +
     *    YHEAD
      CALL PLTVCT(1,XHEAD,YHEAD,A,B)
      A = -ARRLEN*COS(ALPHA)*COS(-THETA) +
     *    ARRLEN*SIN(ALPHA)*SIN(-THETA) + XHEAD
      B = -ARRLEN*COS(ALPHA)*SIN(-THETA) -
     *    ARRLEN*SIN(ALPHA)*COS(-THETA) + YHEAD
      CALL PLTVCT(1,XHEAD,YHEAD,A,B)
      RETURN

      END
