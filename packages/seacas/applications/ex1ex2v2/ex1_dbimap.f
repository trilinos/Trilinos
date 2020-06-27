C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIMAP (NDB, OPTION, NUMEL, MAPEL, *)
C=======================================================================
C$Id: dbimap.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbimap.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:46  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:44  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:10  gdsjaar
c Initial revision
c

C   --*** DBIMAP *** (EXOLIB) Read database element order map
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBIMAP reads the element order map from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'M' to store map
C   --   NUMEL - IN - the number of elements
C   --   MAPEL - OUT - the element order map (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of map upon entry;
C   --upon exit at end of map.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMEL
      INTEGER MAPEL(*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      (MAPEL(IEL), IEL=1,NUMEL)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'ELEMENT ORDER MAP'
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
