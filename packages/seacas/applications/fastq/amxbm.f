C $Id: amxbm.f,v 1.2 1998/07/14 18:18:20 gdsjaar Exp $
C $Log: amxbm.f,v $
C Revision 1.2  1998/07/14 18:18:20  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:03:33  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:31  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]AMXBM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE AMXBM (NPNODE, NPELEM, NXK, AMESUR, BMESUR, KNODE)
C***********************************************************************
C
C  SUBROUTINE AMXBM = ROUTINE TO TRANSFER ELEMENT VARIABLES TO NODES
C
C***********************************************************************
C
      DIMENSION NXK (9, NPELEM)
      DIMENSION AMESUR (NPELEM), BMESUR (NPNODE)
      DIMENSION KNODE (NPNODE)
C
      DO 100 I = 1, NPNODE
         BMESUR(I) = 0.
         KNODE (I) = 0
  100 CONTINUE
C
C  GATHER ALL THE VARIABLES TO THE NODES AND COUNT HOW MANY AT EACH NODE
C
      DO 120 I = 1, NPELEM
         DO 110 J = 1, 4
            NODE = NXK (J, I)
            BMESUR (NODE) = BMESUR (NODE) + AMESUR (I)
            KNODE (NODE) = KNODE (NODE) + 1
  110    CONTINUE
  120 CONTINUE
C
C  GET THE AVERAGE VALUE AT EACH NODE
C
      DO 130 NODE = 1, NPNODE
         BMESUR (NODE) = BMESUR(NODE) / FLOAT (KNODE (NODE))
  130 CONTINUE
C
      RETURN
C
      END
