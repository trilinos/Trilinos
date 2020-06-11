C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOINI (NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL)
C=======================================================================
C$Id: dboini.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dboini.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:24  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:23  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:15  gdsjaar
c Initial revision
c

C   --*** DBOINI *** (EXOLIB) Write database title and initial variables
C   --   Written by Amy Gilkey - revised 12/04/87
C   --
C   --DBOINI writes the title and the initial variables to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   TITLE - IN - the database title
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   NUMEL - IN - the number of elements
C   --   NELBLK - IN - the number of element blocks
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --
C   --Database must be rewound before entry; upon exit positioned at end of
C   --initial variables.

      INTEGER NDB
      CHARACTER*80 TITLE
      INTEGER NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL

      WRITE (NDB) TITLE

      WRITE (NDB) NUMNP, NDIM, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, 1

      RETURN
      END
