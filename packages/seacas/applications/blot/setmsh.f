C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SETMSH (IVIEW, MDEF, MNUM, MLIN, LTYP,
     &   NNPS, ISNPS, NESS, ISESS, MMOD, MTYP,
     &   ISSNPS, ISSESS)
C=======================================================================

C   --*** SETMSH *** (MESH) Set display mode and type
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --SETMSH sets the mesh display options of one or all views.  If the
C   --view is "NONE" or "EMPTY", the other options are set to default
C   --options.
C   --
C   --Parameters:
C   --   IVIEW - IN - the view to be set, 0 for all
C   --   MDEF - IN - the deformed/undeformed option to be set
C   --   MNUM - IN - the view numbering to be set
C   --   MLIN - IN - the display type for the mesh lines to be set
C   --   LTYP - IN - the mesh line type for lines
C   --   NNPS - IN - the number of selected node sets
C   --   ISNPS - IN - the indices of the selected node sets
C   --   NESS - IN - the number of selected side sets
C   --   ISESS - IN - the number of selected side sets
C   --   MMOD - IN - the display mode to be set
C   --   MTYP - IN - the display mode type to be set
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses NUMNPS, NUMESS of /DBNUMG/
C   --   Sets MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET of /MSHOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnumgq.blk'
      include 'mshopt.blk'

      CHARACTER*(*) MDEF, MNUM
      INTEGER LTYP(-1:1)
      INTEGER ISNPS(*)
      INTEGER ISESS(*)
      CHARACTER*(*) MMOD, MTYP
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)

      INTEGER NDEFVW, IXVW

      IF (IVIEW .GE. 1) THEN
         IF ((MDEF .EQ. 'NONE') .OR. (MDEF .EQ. 'EMPTY')) THEN
            MSHDEF(IVIEW) = MDEF
            MSHNUM(IVIEW) = 'NONE'
            MSHLIN(IVIEW) = MSHNON
            CALL INIINT (3, 0, MLNTYP(-1,IVIEW))
            NNPSET(IVIEW) = 0
            CALL INIINT (NUMNPS, 0, ISSNPS(1,IVIEW))
            NESSET(IVIEW) = 0
            CALL INIINT (NUMESS, 0, ISSESS(1,IVIEW))

            CALL SETMOD (IVIEW, 'NONE', ' ')

         ELSE
            MSHDEF(IVIEW) = MDEF
            MSHNUM(IVIEW) = MNUM
            MSHLIN(IVIEW) = MLIN
            CALL CPYINT (3, LTYP(-1), MLNTYP(-1,IVIEW))
            NNPSET(IVIEW) = NNPS
            CALL CPYINT (NNPS, ISNPS, ISSNPS(1,IVIEW))
            NESSET(IVIEW) = NESS
            CALL CPYINT (NESS, ISESS, ISSESS(1,IVIEW))

            IF (MMOD .NE. ' ') CALL SETMOD (IVIEW, MMOD, MTYP)
         END IF

      ELSE
         IF ((MDEF .EQ. 'NONE') .OR. (MDEF .EQ. 'EMPTY')) THEN
            DO 100 I = 1, 4
               MSHDEF(I) = MDEF
               MSHNUM(I) = 'NONE'
               MSHLIN(I) = MSHNON
               CALL INIINT (3, 0, MLNTYP(-1,I))
               NNPSET(I) = 0
               CALL INIINT (NUMNPS, 0, ISSNPS(1,I))
               NESSET(I) = 0
               CALL INIINT (NUMESS, 0, ISSESS(1,I))

               CALL SETMOD (I, 'NONE', ' ')
  100       CONTINUE

         ELSE
            DO 110 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               MSHDEF(I) = MDEF
               MSHNUM(I) = MNUM
               MSHLIN(I) = MLIN
               CALL CPYINT (3, LTYP(-1), MLNTYP(-1,I))
               NNPSET(I) = NNPS
               CALL CPYINT (NNPS, ISNPS, ISSNPS(1,I))
               NESSET(I) = NESS
               CALL CPYINT (NESS, ISESS, ISSESS(1,I))

               IF (MMOD .NE. ' ') CALL SETMOD (I, MMOD, MTYP)
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
