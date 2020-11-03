C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE UNIQUE (LSTSN, NSEG, MAP, ITMP, NUMNIQ, NUMNP)

C***********************************************************************

C     DESCRIPTION:
C       This routine determines the number of unique node numbers in
C       a side set.

C     FORMAL PARAMETERS:
C       LSTSN   INTEGER   List of nodes on this boundary
C       NSEG    INTEGER   Number of nodes in side set
C       MAP     INTEGER   Relates node in side set to list of unique nodes
C       ITMP    INTEGER   Temporary array for sorting nodes
C       NUMNIQ  INTEGER   Number of unique nodes
C       NDIM    INTEGER   Number of spatial dimensions

C     CALLED BY:

C***********************************************************************

      DIMENSION LSTSN(*), MAP(*), ITMP(*)

      CALL INIINT (NSEG, 0, MAP)
      CALL INIINT (NUMNP, 0, ITMP)

      NUMNIQ = 0
      DO 30 I = 1 , NSEG
          IF ( ITMP(LSTSN(I)) .EQ. 0 ) THEN
              NUMNIQ = NUMNIQ + 1
              ITMP(LSTSN(I)) = NUMNIQ
              MAP(I) = NUMNIQ
          ELSE
              MAP(I) = ITMP(LSTSN(I))
          END IF
   30 CONTINUE

      RETURN
      END
