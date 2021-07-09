C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SAVCMD (INLINE, INTYP, CFIELD, NAMES, *)
C=======================================================================

C   --*** SAVCMD *** (ALGEBRA) Perform SAVE command
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --SAVCMD processes the input SAVE command.  It adds all the variables
C   --to be saved to the /VAR../ arrays.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INTYP - IN - the field types
C   --   CFIELD - IN - the character fields
C   --   NAMES - IN - the global, nodal, and element variable names
C   --   * - return statement if command not executed
C   --
C   --Common Variables:
C   --   Sets NUMINP, IXLHS of /VAR../
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_var.blk'
      include 'ag_dbnums.blk'

      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) NAMES(*)

      LOGICAL FFEXST
      CHARACTER*(maxnam) WORD, NAME
      CHARACTER TYPABR
      CHARACTER*5 STRA, STRB

      CHARACTER*(MXSTLN) TYPTBL(5)
      SAVE TYPTBL
      DATA TYPTBL /
     &  'GLOBALS                         ',
     *  'NODALS                          ',
     *  'ELEMENTS                        ',
     *  'ALL                             ',
     &  '                                ' /
c      DATA TYPTBL /
c     &   'HISTORY ', 'GLOBALS ', 'NODALS  ', 'ELEMENTS', 'ALL     ',
c     &   '        ' /

C   --Save the /VAR../ indices so they can be restored in case of error
      NINP = NUMINP
      ILHS = IXLHS

      IF (.NOT. FFEXST (1, INTYP)) THEN
         CALL PRTERR ('CMDERR', 'No options on SAVE command')
         GOTO 120
      END IF

      IFLD = 1
  100 CONTINUE
      IF (FFEXST (IFLD, INTYP)) THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

         IX = LOCSTR (NAME, NVARGL+NVARNP+NVAREL, NAMES)
         IF (IX .GT. 0) THEN
            CALL FFADDC (NAME, INLINE)
            CALL DBVTYP (IX, TYPABR, IVAR)
         ELSE
            IVAR = 0
            if ((name .eq. 'NODE') .or. (name .eq. 'NODES'))
     &         name = 'NODAL'
            CALL ABRSTR (WORD, NAME, TYPTBL)
            TYPABR = WORD(1:1)
            IF (TYPABR .EQ. ' ') THEN
               CALL PRTERR ('CMDWARN', 'Invalid SAVE option "'
     &            // NAME(:LENSTR(NAME)) // '", ignored')
               GOTO 110
            END IF
            CALL FFADDC (WORD, INLINE)
         END IF

         IF ((TYPABR .EQ. 'G') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('G', 1, IGV)
               CALL ADDVAR (NAMES(IGV), NVARGL, 'G', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'G', 1, NINP, ILHS)
            END IF
         END IF

         IF ((TYPABR .EQ. 'N') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('N', 1, INV)
               CALL ADDVAR (NAMES(INV), NVARNP, 'N', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'N', IVAR, NINP, ILHS)
            END IF
         END IF

         IF ((TYPABR .EQ. 'E') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('E', 1, IEV)
               CALL ADDVAR (NAMES(IEV), NVAREL, 'E', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'E', IVAR, NINP, ILHS)
            END IF
         END IF

  110    CONTINUE
         GOTO 100
      END IF

      IF (NUMINP .GE. IXLHS) THEN
         N = NUMINP + (MAXVAR - IXLHS + 1)
         CALL INTSTR (1, 0, N, STRA, LSTRA)
         CALL INTSTR (1, 0, MAXVAR, STRB, LSTRB)
         CALL PRTERR ('CMDSPEC',
     &      'Too many variable names to store, '
     &      // STRA(:LSTRA) // ' > ' // STRB(:LSTRB))
         CALL PRTERR ('CMDERR', 'SAVE command ignored')

         NUMINP = NINP
         IXLHS = ILHS

         GOTO 120
      END IF

      RETURN

  120 CONTINUE
      RETURN 1
      END
