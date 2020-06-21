C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: memext.f,v 1.1 1993/07/16 16:46:58 gdsjaar Exp $
C $Log: memext.f,v $
C Revision 1.1  1993/07/16 16:46:58  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION MEMEXT(IB,MEMRY,INCL)
      INTEGER IB
      INTEGER MEMRY(*)
      INTEGER DELTA
      INTEGER MEMALL
      LOGICAL MEMFRE
      INTEGER LA
      INTEGER IBT
C ... Guess by GDS -- DELTA not defined, set to 0 since that is
C     what it would be on VMS if undefined.
      DELTA = 0
      LA = MEMRY(IB-2) - IB
      IBT = IB
      INCL = LA + DELTA
      IB = MEMALL(INCL,MEMRY)
      DO 2780 I = 1,LA
         MEMRY(IB+I-1) = MEMRY(IBT+I-1)
 2780 CONTINUE
      IF (.NOT.MEMFRE(IBT,MEMRY)) THEN
         WRITE (6,*) 'Error deallocating old segment.'
         MEMEXT = .FALSE.
         RETURN

      END IF

      MEMEXT = .TRUE.
      RETURN

      END
