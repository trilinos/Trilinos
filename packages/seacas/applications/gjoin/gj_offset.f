C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
C $Id: offset.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: offset.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.3  1997/11/18 15:55:30  gdsjaar
C Changes to improve efficiency (factor of 2 to 15 on some Goodyear Treads)
C
C Redid nodeset munching to reduce number of locint calls. First do a
C quick scan on array to get maximum node id in searched array. Then,
C before doing 'locint' call, ensure that id being searched for is less
C than maximum (most times it won't be, saves searching entire list).
C
C Modified the node matching code to also index the nodes within the
C overlap region if doing a nodeset match. This has greatest benefit if
C not all nodes in the nodeset will be matched.
C
C Minor change in offset f -- 'dimension' to 'real'
C
C Revision 1.2  1992/09/03 20:35:02  gdsjaar
C Fixed up handling of mirror and offset - added a scale command which
C is combined with mirror, subroutine offset now does both offset and
C scale.  Added help to irennp.
C
c Revision 1.1  1992/09/02  22:58:05  gdsjaar
c Added mirroring and offsetting capability to gjoin.  Implemented at
c the equivalencing prompt.  Only affects the second database.
c
C=======================================================================
      subroutine offset (off, scale, crd, length)
      real crd(length)

      do 10 i=1, length
         crd(i) = scale * crd(i) + off
 10   continue

      return
      end
