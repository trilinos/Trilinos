C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDMAP (NTXT, NUM, MAP, *)
C=======================================================================

C   --*** RDMAP *** (TXTEXO) Read database element order map
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDMAP reads the element order map, node number map, and element
C     number maps from the text file.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NUM - IN - the number of elements or nodes
C   --   MAP - OUT - the map
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of map upon entry;
C   --upon exit at end of map.

      INTEGER MAP(*)
      CHARACTER*80 LINE

      READ (NTXT, *, END=100, ERR=100)
      READ (NTXT, '(A)', END=100, ERR=100) LINE
C ... Check for keyword 'sequence' at beginning of line
      if (line(:8) .eq. 'sequence' .or. line(:8) .eq. 'SEQUENCE') then
        do 10 i=1, num
          map(i) = i
 10     continue
      else
        READ (NTXT, *, END=100, ERR=100) (MAP(I), I=1,NUM)
      end if

      RETURN

  100 CONTINUE
      CALL PRTERR ('FATAL', 'Reading MAP')
      RETURN 1
      END
