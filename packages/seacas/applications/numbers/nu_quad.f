C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QUAD(XXX, XI, XG, NDIM, NNODES, NQUAD, WT)
      DIMENSION XXX(NDIM+1,NNODES,NQUAD), XI(NDIM,*), XG(NDIM,*)

      IF (NQUAD .EQ. 1) THEN
          QUADL = 0.0
      ELSE
          QUADL = 1./SQRT(3.)
      END IF

      WT = 2.**NDIM / DBLE(NQUAD)
      IF (NQUAD .EQ. 1) THEN
          XG(1,1) = 0.0
          XG(2,1) = 0.0
          XG(3,1) = 0.0
      ELSE
          DO 20 I=1, NNODES
              DO 10 J=1, NDIM
                  XG(J,I) = XI(J,I) * QUADL
   10         CONTINUE
   20     CONTINUE
      END IF

      IF (NDIM .EQ. 3) THEN
          DO 40 I=1, NQUAD
              DO 30 J=1, NNODES
                  TMP1 = (1. + XI(1,J) * XG(1,I))
                  TMP2 = (1. + XI(2,J) * XG(2,I))
                  TMP3 = (1. + XI(3,J) * XG(3,I))

                  XXX(1,J,I) = TMP1    * TMP2 * TMP3 / 8.0
                  XXX(2,J,I) = XI(1,J) * TMP2 * TMP3 / 8.0
                  XXX(3,J,I) = XI(2,J) * TMP1 * TMP3 / 8.0
                  XXX(4,J,I) = XI(3,J) * TMP1 * TMP2 / 8.0
   30         CONTINUE
   40     CONTINUE
      ELSE
          DO 60 I=1, NQUAD
              DO 50 J=1, NNODES
                  TMP1 = (1. + XI(1,J) * XG(1,I))
                  TMP2 = (1. + XI(2,J) * XG(2,I))

                  XXX(1,J,I) = TMP1    * TMP2 / 4.0
                  XXX(2,J,I) = XI(1,J) * TMP2 / 4.0
                  XXX(3,J,I) = XI(2,J) * TMP1 / 4.0
   50         CONTINUE
   60     CONTINUE
      END IF
      RETURN
      END
