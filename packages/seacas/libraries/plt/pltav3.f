C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltav3.f,v 1.4 2000/10/25 18:55:02 gdsjaar Exp $
C $Log: pltav3.f,v $
C Revision 1.4  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.3  1993/07/19 17:06:32  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:09  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:41  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTAV3(UMAP,N,X1,Y1,Z1,X2,Y2,Z2,TH,XL)
      REAL UMAP(*)
      INTEGER N
      REAL X1(*),Y1(*),Z1(*)
      REAL X2(*),Y2(*),Z2(*)
      REAL TH
      REAL XL
      REAL PX(32),PY(32),QX(32),QY(32)
      INTEGER MASK(1)
      include 'izbit.inc'

      MASK(1) = -1
      CALL PLTMV3(UMAP,N,MASK,X1,Y1,Z1,X2,Y2,Z2,PX,PY,QX,QY)
      DO 2020 I = 1,N
         IF (IAND(MASK(1),IZBIT(I)).NE.0) THEN
            CALL PLTARR(PX(I),PY(I),QX(I),QY(I),TH,XL)
         END IF

 2020 CONTINUE
      RETURN

      END
