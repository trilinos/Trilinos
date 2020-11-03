C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETCAV (ERROR, IDESS, NUMESS)
      DIMENSION IDESS(*)
      include 'nu_cav.blk'
      LOGICAL ERROR
      CHARACTER*80 STRA

      PARAMETER (MXFLD = 12)
      DIMENSION RV(MXFLD), KV(MXFLD)
      CHARACTER*32 CV(MXFLD), PRMPT
      PRMPT = '  Cavity Side Set > '
      MAXF = MIN (MXFLD, MAXCAV)
      IF (ICAV(1) .EQ. 0 .OR. NUMCAV .EQ. 0) THEN
         CALL FREFLD (0, 0, PRMPT(:LENSTR(PRMPT)+1), MAXF, IOS,
     *      NUMCAV, KV, CV, ICAV, RV)
      END IF

      ERROR = .FALSE.
      DO 10 NCAV = 1, NUMCAV
         IFND(NCAV) = LOCINT (ICAV(NCAV), NUMESS, IDESS)
         IF (IFND(NCAV) .EQ. 0) THEN
            WRITE (STRA, 20) ICAV(NCAV)
            CALL SQZSTR (STRA, LSTR)
            CALL PRTERR ('ERROR', STRA(:LSTR))
            ERROR = .TRUE.
         END IF
   10 CONTINUE
   20 FORMAT (' Cavity Flag ',I5,' not found. ')

      RETURN
      END
