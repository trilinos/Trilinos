C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE AMXBM (NPNODE, NPELEM, NXK, AMESUR, BMESUR, KNODE)
C***********************************************************************

C  SUBROUTINE AMXBM = ROUTINE TO TRANSFER ELEMENT VARIABLES TO NODES

C***********************************************************************

      DIMENSION NXK (9, NPELEM)
      DIMENSION AMESUR (NPELEM), BMESUR (NPNODE)
      DIMENSION KNODE (NPNODE)

      DO 100 I = 1, NPNODE
         BMESUR(I) = 0.
         KNODE (I) = 0
  100 CONTINUE

C  GATHER ALL THE VARIABLES TO THE NODES AND COUNT HOW MANY AT EACH NODE

      DO 120 I = 1, NPELEM
         DO 110 J = 1, 4
            NODE = NXK (J, I)
            BMESUR (NODE) = BMESUR (NODE) + AMESUR (I)
            KNODE (NODE) = KNODE (NODE) + 1
  110    CONTINUE
  120 CONTINUE

C  GET THE AVERAGE VALUE AT EACH NODE

      DO 130 NODE = 1, NPNODE
         BMESUR (NODE) = BMESUR(NODE) / DBLE(KNODE (NODE))
  130 CONTINUE

      RETURN

      END
