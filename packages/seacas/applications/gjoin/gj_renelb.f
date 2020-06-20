C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENELB (NELBLK, IOFFNP, IXNP, NUMELB, NUMLNK, LINK)
C=======================================================================
C $Id: renelb.f,v 1.1 1999/01/18 19:21:26 gdsjaar Exp $
C $Log: renelb.f,v $
C Revision 1.1  1999/01/18 19:21:26  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:00  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:58  gdsjaar
c Initial revision
c

C   --*** RENELB *** (GJOIN) Renumber connectivity in element blocks
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENELB renumbers the nodes in the connectivity arrays.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   IOFFNP - IN - the nodal offset: if positive, add to node number;
C   --      if negative, use IXNP
C   --   IXNP - IN - the new node number for each node
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK - IN - the connectivity for each block

      INTEGER IXNP(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)

      ILNK = 1

      DO 100 IELB = 1, NELBLK
         CALL RENIX (NUMELB(IELB) * NUMLNK(IELB),
     &      IOFFNP, IXNP, LINK(ILNK))
         ILNK = ILNK + NUMLNK(IELB) * NUMELB(IELB)
  100 CONTINUE

      RETURN
      END
