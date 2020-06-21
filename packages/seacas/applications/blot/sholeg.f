C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: sholeg.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:12:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/09/24  17:32:43  gdsjaar
c Added an outline off/on command to toggle drawing of the view window outline
c
c Revision 1.2  1990/12/14  08:57:22  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SHOLEG (SHOTYP, XYTYPE, DOQA, DOLEG, DOAXIS, CAPTN,
     *  DOBOX)
C=======================================================================

C   --*** SHOLEG *** (BLOT) Display plot labeling option
C   --   Written by Amy Gilkey - revised 04/28/88
C   --
C   --SHOLEG display the plot labeling option requested by the option type:
C   --   QA - whether the QA information is to be drawn on legend
C   --   LEGEND - whether the non-QA legend information is to be drawn
C   --   AXIS - whether the axis is to be numbered
C   --   CAPTION - the plot caption
C   --
C   --Parameters:
C   --   SHOTYP - IN - the show option (see above)
C   --   XYTYPE - IN - true iff current program is an XY curve versus mesh plot
C   --   DOQA - IN - true iff QA information is to be drawn on legend
C   --   DOLEG - IN - true iff non-QA legend information is to be drawn
C   --   DOAXIS - IN - true iff axis is to be labeled
C   --   CAPTN - IN - the three-line plot caption
C   --   DOBOX - IN - true iff plot is to be outlined

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) SHOTYP
      LOGICAL XYTYPE
      LOGICAL DOQA, DOLEG, DOAXIS, DOBOX
      CHARACTER*(*) CAPTN(3)

      CHARACTER*20 STR20

      IF (.NOT. XYTYPE) THEN
         STR20 = 'mesh'
      ELSE
         STR20 = 'X-Y curve'
      END IF
      LSTR = LENSTR (STR20)

      IF (SHOTYP .EQ. 'QA') THEN
         IF (DOQA) THEN
            WRITE (*, 10000) 'Include QA information on ',
     &         STR20(:LSTR), ' plot legend'
         ELSE
            WRITE (*, 10000) 'Omit QA information from ',
     &         STR20(:LSTR), ' plot legend'
         END IF

      ELSE IF (SHOTYP .EQ. 'LEGEND') THEN
         IF (DOLEG) THEN
            WRITE (*, 10000) 'Include non-QA legend information on ',
     &         STR20(:LSTR), ' plot'
         ELSE
            WRITE (*, 10000) 'Omit non-QA legend information on ',
     &         STR20(:LSTR), ' plot'
         END IF

      ELSE IF (SHOTYP .EQ. 'AXIS') THEN
         IF (DOAXIS) THEN
            WRITE (*, 10000) 'Number ', STR20(:LSTR), ' axes'
         ELSE
            WRITE (*, 10000) 'Do not number ', STR20(:LSTR), ' axes'
         END IF

      ELSE IF (SHOTYP .EQ. 'OUTLINE') THEN
         IF (DOBOX) THEN
            WRITE (*, 10000) 'Outline Plot'
         ELSE
            WRITE (*, 10000) 'Do not outline plot'
         END IF

      ELSE IF (SHOTYP .EQ. 'CAPTION') THEN
         DO 100 IEND = 3, 1, -1
            IF (CAPTN(IEND) .NE. ' ') GOTO 110
  100    CONTINUE
  110    CONTINUE
         IF (STR20 .EQ. 'mesh') STR20 = 'Mesh'
         IF (IEND .LE. 0) THEN
            WRITE (*, 10000)
     &         STR20(:LSTR), ' plot caption is not defined'
         ELSE
            WRITE (*, 10000) STR20(:LSTR), ' plot caption:'
            DO 120 I = 1, IEND
               WRITE (*, 10000) '   ', CAPTN(I)(:LENSTR(CAPTN(I)))
  120       CONTINUE
         END IF
      END IF

      RETURN
10000  FORMAT (1X, 5A)
      END
