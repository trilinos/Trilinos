C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: chrlen.f,v $
C Revision 1.3  1993/07/19 14:18:12  gdsjaar
C Reformatted flow of control
C
c Revision 1.2  1993/07/19  14:16:17  gdsjaar
c Reformatted flow of control
c
c Revision 1.1  1993/07/16  16:46:20  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      INTEGER FUNCTION CHRLEN(S)
      CHARACTER*(*) S

      L = LEN(S)
      J = INDEX(S,CHAR(0))
      IF (J.EQ.0) THEN
         I = L

      ELSE

         I = J - 1
      END IF

   10 CONTINUE

      IF ((I.GT.0.AND.S(I:I).EQ.' ')) THEN
         I = I - 1
         GO TO 10

      END IF

      CHRLEN = I

      END
