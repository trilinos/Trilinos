C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE CMDMOD (VERB, VERB2, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   IVIEW, JVIEW, NAMES, NVOLD, NCOLD, FIXCON, VECSCL,
     &   ISSNPS, ISSESS, *)
C=======================================================================

C   --*** CMDMOD *** (DETOUR) Process display mode commands
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --Parameters:
C   --   VERB, VERB2 - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   IVIEW - IN - the view number
C   --   JVIEW - IN - IVIEW (if non-zero) or a defined non-empty view number
C   --      (if any)
C   --   NAMES - IN - the variable names
C   --   NVOLD - IN/OUT - the old selected variables
C   --   NCOLD - IN/OUT - the old contour variable
C   --   FIXCON - IN/OUT - true iff contour parameters are fixed
C   --   VECSCL - IN/OUT - the vector scale factor
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses NVARNP, NVAREL of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Sets and uses MSHDEF, MSHNUM, MSHLIN, MLNTYP of /MSHOPT/
C   --   Sets and uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'mshopt.blk'
      include 'detopt.blk'

      CHARACTER*(*) VERB, VERB2
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMES(*)
      INTEGER NVOLD(4), NCOLD
      LOGICAL FIXCON
      REAL VECSCL
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      LOGICAL FFEXST, MATSTR
      INTEGER NUMMOD
      CHARACTER*(MXNAME) WORD
      CHARACTER*(MXNAME) MMOD, MTYP
      INTEGER LTYP(-1:1)
      LOGICAL LDUM

C   --Get the appropriate verb for a variable, based on the display mode

      MTYP = ' '
      IF (VERB .EQ. ' ') THEN
         IF (MODDET(JVIEW) .EQ. 'CONTOUR') THEN
            IF (MODTYP(JVIEW) .EQ. 'LINE') THEN
               VERB = 'CONTOUR'
            ELSE
               VERB = 'PAINT'
            END IF
         ELSE IF (MODDET(JVIEW) .EQ. 'CONTOUR') THEN
            IF (MODTYP(JVIEW) .EQ. 'PAINT') THEN
               VERB = 'EPAINT'
            END IF
         ELSE IF (MODDET(JVIEW) .EQ. 'VECTOR') THEN
            VERB = 'VECTOR'
         ELSE IF (MODDET(JVIEW) .EQ. 'SYMBOL') THEN
            VERB = 'SYMBOL'
            MTYP = MODTYP(JVIEW)
         ELSE IF (MODDET(JVIEW) .EQ. 'GAUSS') THEN
            VERB = 'GAUSS'
            MTYP = MODTYP(JVIEW)
         ELSE
            VERB = 'CONTOUR'
         END IF
      END IF

      CALL INIINT (3, 1, LTYP)

      IF (VERB .EQ. 'WIREFRAM') THEN
         CALL FFADDC (VERB, INLINE)
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHSEL, LTYP,
     &      0, IDUM, 0, IDUM, 'WIREFRAM', ' ', ISSNPS, ISSESS)

      ELSE IF (VERB .EQ. 'SOLID') THEN
         CALL FFADDC (VERB, INLINE)
         LTYP(1) = -1
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHSEL, LTYP,
     &      0, IDUM, 0, IDUM, 'SOLID', ' ', ISSNPS, ISSESS)

      ELSE IF ((VERB .EQ. 'CONTOUR') .OR. (VERB .EQ. 'PAINT')) THEN
         CALL FFADDC (VERB, INLINE)
         IF (VERB .EQ. 'CONTOUR') THEN
            MTYP = 'LINE'
         ELSE
            MTYP = 'PAINT'
         END IF

C      --If variable given, force the contour range to be re-calculated and
C      --set the contour flag to undetermined
         IF (FFEXST (IFLD, INTYP)) THEN
            NCOLD = 0
            FIXCON = .FALSE.
         END IF

         CALL CMDVAR ('CONTOUR', MTYP, NAMES,
     &      INLINE, IFLD, INTYP, CFIELD, IDTVAR, *100)

         LTYP(1) = 2
         IF (MTYP .EQ. 'PAINT') LTYP(1) = -LTYP(1)
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHDIV, LTYP,
     &      0, IDUM, 0, IDUM, 'CONTOUR', MTYP, ISSNPS, ISSESS)

         IF ((.NOT. FIXCON) .AND. (IVIEW .NE. 0)) THEN
            IF (((NUMMOD (MODDET, MODTYP, 'CONTOUR', 'LINE')
     &         + NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'LINE')) .GE. 1)
     &         .AND. ((NUMMOD (MODDET, MODTYP, 'CONTOUR', 'PAINT')
     &         + NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'PAINT')) .GE. 1))
     &         THEN
               CALL PRTERR ('CMDWARN', 'Contradictory contour types')
            END IF
         END IF

         CALL CHKVAR (MODDET, MODTYP, IVIEW, IDTVAR, NVOLD,
     &      NNDVAR, NEDVAR, LDUM)

      ELSE IF (VERB .EQ. 'EPAINT') THEN
         CALL FFADDC (VERB, INLINE)
         IF (VERB .EQ. 'EPAINT') THEN
            MTYP = 'PAINT'
         END IF

C      --If variable given, force the contour range to be re-calculated and
C      --set the contour flag to undetermined
         IF (FFEXST (IFLD, INTYP)) THEN
            NCOLD = 0
            FIXCON = .FALSE.
         END IF

         CALL CMDVAR ('ELEMCONT', MTYP, NAMES,
     &      INLINE, IFLD, INTYP, CFIELD, IDTVAR, *100)

         LTYP(1) = 2
         IF (MTYP .EQ. 'PAINT') LTYP(1) = -LTYP(1)
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHDIV, LTYP,
     &      0, IDUM, 0, IDUM, 'ELEMCONT', MTYP, ISSNPS, ISSESS)

         IF ((.NOT. FIXCON) .AND. (IVIEW .NE. 0)) THEN
            IF (((NUMMOD (MODDET, MODTYP, 'CONTOUR', 'LINE')
     &         + NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'LINE')) .GE. 1)
     &         .AND. ((NUMMOD (MODDET, MODTYP, 'CONTOUR', 'PAINT')
     &         + NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'PAINT')) .GE. 1))
     &         THEN
               CALL PRTERR ('CMDWARN', 'Contradictory contour types')
            END IF
         END IF

         CALL CHKVAR (MODDET, MODTYP, IVIEW, IDTVAR, NVOLD,
     &      NNDVAR, NEDVAR, LDUM)

      ELSE IF (VERB .EQ. 'VECTOR') THEN
         CALL FFADDC (VERB, INLINE)
C      --The vector type (nodal or element) is set in CHKVAR
         MTYP = ' '

         CALL CMDVAR ('VECTOR', MTYP, NAMES,
     &      INLINE, IFLD, INTYP, CFIELD, IDTVAR, *100)

         LTYP(1) = 2
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHDIV, LTYP,
     &      0, IDUM, 0, IDUM, 'VECTOR', MTYP, ISSNPS, ISSESS)

C      --Reset the vector/symbol scaling factor
         VECSCL = 1.0
         VERB2 = 'VECSCL'

         CALL CHKVAR (MODDET, MODTYP, IVIEW, IDTVAR, NVOLD,
     &      NNDVAR, NEDVAR, LDUM)

      ELSE IF ((VERB .EQ. 'SYMBOL') .OR. (VERB .EQ. 'GAUSS')) THEN
         CALL FFADDC (VERB, INLINE)
         IF (VERB .EQ. 'SYMBOL') THEN
            MMOD = 'SYMBOL'
         ELSE
            MMOD = 'GAUSS'
         END IF

         IF (MTYP .EQ. ' ') THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ANGLE', 1)) THEN
               CALL FFADDC ('ANGLE', INLINE)
               MTYP = 'ANGLE'
            ELSE IF (MATSTR (WORD, 'CRACK', 1)) THEN
               CALL FFADDC ('CRACK', INLINE)
               MTYP = 'CRACK'
            ELSE IF (MATSTR (WORD, 'STATE', 2)) THEN
               CALL FFADDC ('STATE', INLINE)
               MTYP = 'STATE'
            ELSE IF (MATSTR (WORD, 'USER', 1)) THEN
               CALL FFADDC ('USER', INLINE)
               MTYP = 'USER'
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "ANGLE", "CRACK", OR "STATE"'
     &            // ' before variable name')
               GOTO 100
            END IF
         END IF

         CALL CMDVAR (MMOD, MTYP, NAMES,
     &      INLINE, IFLD, INTYP, CFIELD, IDTVAR, *100)

         LTYP(1) = 2
         CALL SETMSH (IVIEW, 'DEFORM', 'NONE', MSHDIV, LTYP,
     &      0, IDUM, 0, IDUM, MMOD, MTYP, ISSNPS, ISSESS)

         IF (MTYP .NE. 'SPHERE'  .AND.  MTYP .NE. 'FSPHERE') THEN
C           --Reset the vector/symbol scaling factor
            VECSCL = 1.0
            VERB2 = 'VECSCL'
            CALL CHKVAR (MODDET, MODTYP, IVIEW, IDTVAR, NVOLD,
     &         NNDVAR, NEDVAR, LDUM)
         ENDIF

      ELSE
         GOTO 100
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
