C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: cmdvwc.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:56:31  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:37  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CMDVWC (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   IVIEW, JVIEW, *)
C=======================================================================

C   --*** CMDVWC *** (MESH) Process VIEW command (header only)
C   --   Written by Amy Gilkey - revised 04/13/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the VIEW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   IVIEW - IN - the view number
C   --   JVIEW - IN - IVIEW (if non-zero) or a defined non-empty view number
C   --      (if any)
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'mshopt.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)

      INTEGER NDEFVW, IXVW

      IF (VERB .EQ. 'VIEW') THEN
         CALL FFADDC (VERB, INLINE)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'view number', -1, IVIEW, *100)
         CALL FFADDI (IVIEW, INLINE)
         IF ((IVIEW .LT. 0) .OR. (IVIEW .GT. 4)) THEN
            CALL PRTERR ('CMDERR', 'Expected view number')
            GOTO 100
         END IF
         IF ((IVIEW .GE. 1) .AND. (MSHDEF(IVIEW) .EQ. 'NONE')) THEN
            CALL PRTERR ('CMDERR', 'Specified view is undefined')
            GOTO 100
         END IF
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', VERB)

      ELSE
         IF (NDEFVW (.TRUE.) .GT. 1) THEN
            CALL PRTERR ('CMDERR',
     &         'VIEW command must be used with multiple views')
            GOTO 100
         END IF

         IVIEW = 2
      END IF

      IF (IVIEW .EQ. 0) THEN
         JVIEW = IXVW (.FALSE., 1)
         IF (JVIEW .EQ. 0) THEN
            CALL PRTERR ('CMDERR',
     &         'VIEW 0 not allowed because all views are empty')
            GOTO 100
         END IF
      ELSE
         JVIEW = IVIEW
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
