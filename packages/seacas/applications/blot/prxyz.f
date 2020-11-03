C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRXYZ (OPTION, NOUT, NDIM, NAMECO, NUMNP, LISNP,
     &   XN, YN, ZN, MAPND)
C=======================================================================

C   --*** PRXYZ *** (BLOT) Display database coordinates
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --PRXYZ displays the coordinate array.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NDIM - IN - the number of coordinates per node
C   --   NAMECO - IN - the coordinate names
C   --   NUMNP - IN - the number of nodes
C   --   LISNP - IN - the indices of the selected nodes
C   --   XN, YN, ZN - IN - the nodal coordinates

      CHARACTER*(*) OPTION
      CHARACTER*(*) NAMECO(*)
      INTEGER LISNP(0:*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER MAPND(*)

      LOGICAL ISABRT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) (NAMECO(I)(:8), I=1,NDIM)
      ELSE
         WRITE (*, 10010) (NAMECO(I)(:8), I=1,NDIM)
      END IF

      DO 100 IX = 1, LISNP(0)
         IF (ISABRT ()) RETURN
         INP = LISNP(IX)
         IF (NOUT .GT. 0) THEN
            IF (NDIM .LE. 2) THEN
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &          MAPND(INP), XN(INP), YN(INP)
            ELSE
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &          MAPND(INP), XN(INP), YN(INP), ZN(INP)
            END IF
         ELSE
            IF (NDIM .LE. 2) THEN
               WRITE (*, 10020, IOSTAT=IDUM)
     &          MAPND(INP), XN(INP), YN(INP)
            ELSE
               WRITE (*, 10020, IOSTAT=IDUM)
     &          MAPND(INP), XN(INP), YN(INP), ZN(INP)
            END IF
         END IF
  100 CONTINUE

      RETURN

10000  FORMAT (/, 1X, 'COORDINATES')
10010  FORMAT (/, 1X, 'Node Ids are Global Ids',/,
     *   5X, 6X, 4X, 5 (3X, A8, :, 3X), :, /,
     &   (1X, 5 (3X, A8, :, 3X)))
10020  FORMAT (1X, 'Node ', I9, 4X, 5 (1X, 1PE13.6), :, /,
     &   (1X, 5 (1X, 1PE13.6)))
      END
