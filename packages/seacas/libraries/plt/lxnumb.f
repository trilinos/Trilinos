C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxnumb.f,v 1.2 1993/07/16 18:28:50 gdsjaar Exp $
C $Log: lxnumb.f,v $
C Revision 1.2  1993/07/16 18:28:50  gdsjaar
C Changed real*8 to double precision
C
c Revision 1.1  1993/07/16  16:46:44  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      LOGICAL FUNCTION LXNUMB(N,ND,CH)
      CHARACTER*(*) CH
      DOUBLE PRECISION N
      LOGICAL LXSET

      N = 0.
      ND = 0
 2580 IF (.NOT. (LXSET('0123456789',CH))) GO TO 2590
      N = N*10. + DBLE(ICHAR(CH)-ICHAR('0'))
      ND = ND + 1
      GO TO 2580

 2590 CONTINUE
      IF (ND.EQ.0) THEN
         LXNUMB = .FALSE.
         RETURN

      END IF

      LXNUMB = .TRUE.
      RETURN

      END
