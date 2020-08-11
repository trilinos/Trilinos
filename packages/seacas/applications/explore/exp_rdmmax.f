C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDMMAX (IFLD, INTYP, CFIELD,
     &     NAMIGV, NAMINV, NAMIEV,
     &     NAMOGV, NAMONV, NAMOEV,
     &     NSTEP, MMSTEP, MMNAME, MMTYP, MMVAR, MMNUM, *)
C=======================================================================

C     --*** RDMMAX *** (EXPLORE) Parse the MINMAX command parameters
C     --
C     --RDMMAX parses the MINMAX command parameters.  It first reads
C     --the name, if any, and determines its type and variable number.
C     --It then reads the 'ALL' versus 'THIS' parameter.
C     --
C     --Parameters:
C     --   IFLD - IN/OUT - the number of the next entry to scan, incremented
C     --   INTYP - IN - the free-format field types
C     --      -1 = none, 0 = name, 1 = real, 2 = integer
C     --   CFIELD - IN - the input character fields
C     --   NAMIGV - IN - the names of the global variables as input
C     --   NAMINV - IN - the names of the nodal variables as input
C     --   NAMIEV - IN - the names of the element variables as input
C     --   NAMOGV - IN - the names of the global variables for comparison
C     --   NAMONV - IN - the names of the nodal variables for comparison
C     --   NAMOEV - IN - the names of the element variables for comparison
C     --   NSTEP - IN - the current step number
C     --   MMSTEP - IN/OUT - the requested step number, <=0 for all
C     --   MMNAME - IN/OUT - min/max variable name
C     --   MMTYP - IN/OUT - min/max variable type:
C     --      'G'lobal, 'N'odal, 'E'lement
C     --   MMVAR - IN/OUT - min/max variable number
C     --   MMNUM - IN/OUT - number of sequential min/max requests for this
C     --      variable
C     --   * - return statement if invalid range of integers specified;
C     --      message printed
C     --
C     --Common Variables:
C     --   Uses NVARNP, NVAREL, NVARGL of /DBNUMS/

      include 'exodusII.inc'
      include 'exp_dbnums.blk'

      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) NAMIGV(*), NAMINV(*), NAMIEV(*)
      CHARACTER*(*) NAMOGV(*), NAMONV(*), NAMOEV(*)
      CHARACTER*(*) MMNAME
      CHARACTER MMTYP

      CHARACTER*(256) NAME, WORD
      CHARACTER CH

      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

      IF (NAME .NE. ' ') THEN

C     --Get the variable name

         IV = LOCSTR (NAME, NVARGL, NAMOGV)
         IF (IV .GT. 0) THEN
            CH = 'G'
            NAME = NAMIGV(IV)
         ELSE
            IV = LOCSTR (NAME, NVARNP, NAMONV)
            IF (IV .GT. 0) THEN
               CH = 'N'
               NAME = NAMINV(IV)
            ELSE
               IV = LOCSTR (NAME, NVAREL, NAMOEV)
               IF (IV .GT. 0) THEN
                  CH = 'E'
                  NAME = NAMIEV(IV)
               END IF
            END IF
         END IF
         IF (IV .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 100
         END IF

C     --Get the ALL (versus this time step) field

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'THIS', WORD)
         IF (WORD .EQ. 'ALL') THEN
            MMSTEP = 0
         ELSE IF (WORD .EQ. 'THIS') THEN
            MMSTEP = NSTEP
         ELSE
            CALL PRTERR ('CMDERR', 'Expected "ALL" or "THIS"')
            GOTO 100
         END IF

         MMNAME = NAME
         MMVAR = IV
         MMTYP = CH
         MMNUM = 1
      ELSE

C     --Get next min/max for last variable

         IF (MMVAR .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 100
         END IF
      END IF

      RETURN

 100  CONTINUE
      RETURN 1
      END
