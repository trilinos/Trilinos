C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOMAP (NDB, NUMEL, MAPEL)
C=======================================================================
C$Id: dbomap.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbomap.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:27  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:25  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:15  gdsjaar
c Initial revision
c

C   --*** DBOMAP *** (EXOLIB) Write database element order map
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBOMAP writes the element order map to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMEL - IN - the number of elements
C   --   MAPEL - IN - the element order map
C   --
C   --Database must be positioned at start of map upon entry;
C   --upon exit at end of map.

      INTEGER NDB
      INTEGER NUMEL
      INTEGER MAPEL(*)

      IF (NUMEL .GT. 0) THEN
         WRITE (NDB) (MAPEL(IEL), IEL=1,NUMEL)
      ELSE
         WRITE (NDB) 0
      END IF

      RETURN
      END
