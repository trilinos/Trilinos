C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE SELMAP (NOUT, TYPE, IDGLO, NUMEL, MAPEL)
C=======================================================================

C   --SELMAP locates the local node corresponding to the global id
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   TYPE - IN - 'node' or 'element'
C   --   NUMEL - IN - the number of node/elements
C   --   MAPEL - IN - the node/element number map

C ... TYPE = Node or Element
      CHARACTER*(*) TYPE
      INTEGER MAPEL(*)

      CHARACTER*128 STRA

      LOGICAL FOUND

      found = .FALSE.
      do 80 i = 1, numel
        if (idglo .eq. mapel(i)) then
          found = .TRUE.
          IF (NOUT .GT. 0) THEN
            write (nout, 10030) TYPE, IDGLO, TYPE, I
          ELSE
            write (*, 10030) TYPE, IDGLO, TYPE, I
          END IF
          go to 90
        end if
 80   continue
 90   continue
      if (.not. found) then
        write (stra, 10000) type, idglo
        call sqzstr(stra, lstra)
        call prterr('WARNING', stra(:lstra))
      end if
      RETURN

10000  FORMAT (' No local ',A,' has global id equal to ',I12)
10030  format (1x, 3x, 'Global ',A,I12,' is local ',A,I12)
      END

