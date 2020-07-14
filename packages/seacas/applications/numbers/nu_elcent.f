C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: elcent.f,v 1.1 1991/02/21 15:43:04 gdsjaar Exp $
C $Log: elcent.f,v $
C Revision 1.1  1991/02/21 15:43:04  gdsjaar
C Initial revision
C
      SUBROUTINE ELCENT (ELCEN, IX, COORD, NDIM, NUMEL, NELNOD, NUMNP)
      DIMENSION ELCEN(NUMEL, *), IX(NELNOD, *), COORD(NUMNP, *)
C
      IF (NDIM .EQ. 2) THEN
          DO 10 I=1, NUMEL
              ELCEN(I,1) = (COORD(IX(1,I),1) + COORD(IX(2,I),1) +
     *                      COORD(IX(3,I),1) + COORD(IX(4,I),1))/4.0
              ELCEN(I,2) = (COORD(IX(1,I),2) + COORD(IX(2,I),2) +
     *                      COORD(IX(3,I),2) + COORD(IX(4,I),2))/4.0
   10     CONTINUE
      ELSE
          DO 20 I=1, NUMEL
              ELCEN(I,1) = (COORD(IX(1,I),1) + COORD(IX(2,I),1) +
     *                      COORD(IX(3,I),1) + COORD(IX(4,I),1) +
     *                      COORD(IX(5,I),1) + COORD(IX(6,I),1) +
     *                      COORD(IX(7,I),1) + COORD(IX(8,I),1))/8.0
C
              ELCEN(I,2) = (COORD(IX(1,I),2) + COORD(IX(2,I),2) +
     *                      COORD(IX(3,I),2) + COORD(IX(4,I),2) +
     *                      COORD(IX(5,I),2) + COORD(IX(6,I),2) +
     *                      COORD(IX(7,I),2) + COORD(IX(8,I),2))/8.0
C
              ELCEN(I,3) = (COORD(IX(1,I),3) + COORD(IX(2,I),3) +
     *                      COORD(IX(3,I),3) + COORD(IX(4,I),3) +
     *                      COORD(IX(5,I),3) + COORD(IX(6,I),3) +
     *                      COORD(IX(7,I),3) + COORD(IX(8,I),3))/8.0
   20     CONTINUE
      END IF
      RETURN
      END
