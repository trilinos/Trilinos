C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRMAP (OPTION, NOUT, TYPE, NUMEL, MAPEL)
C=======================================================================

C   --*** PRMAP *** (EXPLORE) Display database node/element number map
C   --
C   --PRMAP displays the node/element number map.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMEL - IN - the number of node/elements
C   --   MAPEL - IN - the node/element number map

      CHARACTER*(*) OPTION

C ... TYPE = Node or Element
      CHARACTER*(*) TYPE
      INTEGER MAPEL(*)

      CHARACTER*4 FMT
      CHARACTER*32 STRA, STRB

C ... True if map is one-to-one
      LOGICAL MAPONE

      IF (NOUT .GT. 0) WRITE (NOUT, 10000) TYPE

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) TYPE
      ELSE
         WRITE (*, 10010) TYPE
      END IF

      WRITE (STRA, '(I12)', IOSTAT=IDUM) NUMEL
      CALL SQZSTR (STRA, LSTRA)
      WRITE (FMT, '(''(I'', I1, '')'')', IOSTAT=IDUM) LSTRA

C ... Check for 1-1 mapping
      mapone = .TRUE.
      do 80 i = 1, numel
        if (mapel(i) .ne. i) then
          mapone = .FALSE.
          go to 90
        end if
 80   continue
 90   continue
      if (mapone) then
        IF (NOUT .GT. 0) THEN
          write (nout, 10030) TYPE
        ELSE
          write (*, 10030) TYPE
        END IF
      else
        DO 100 IEL = 1, NUMEL, 8
          NE = MIN (IEL+7, NUMEL)
          WRITE (STRA, FMT, IOSTAT=IDUM) IEL
          WRITE (STRB, FMT, IOSTAT=IDUM) NE
          IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020, IOSTAT=IDUM)
     &        STRA(:LSTRA), STRB(:LSTRA), (MAPEL(N), N=IEL,NE)
          ELSE
            WRITE (*, 10020, IOSTAT=IDUM)
     &        STRA(:LSTRA), STRB(:LSTRA), (MAPEL(N), N=IEL,NE)
          END IF
 100    CONTINUE
      end if

      RETURN

10000  FORMAT (/, 1X, A,' NUMBER MAP')
10010  FORMAT (/, 1X, A,' Number Map:')
10020  FORMAT (1X, 3X, A, '..', A, 3X, 8I12)
10030  format (1x, 3x, 'Map does not modify local ',
     &   A, ' ids (X maps to X)')
      END
