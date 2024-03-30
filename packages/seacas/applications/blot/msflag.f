C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSFLAG (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL,
     &   MINMSH, MAXMSH, MAXHID)
C=======================================================================

C   --*** MSFLAG *** (MESH) Set mesh plot flags
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --MSFLAG sets the flags needed to plot a mesh.
C   --
C   --Parameters:
C   --   ANYDEF - OUT - true iff any deformed mesh is to be plotted
C   --   ANYUND - OUT - true iff any undeformed mesh is to be plotted
C   --   DOIXF - OUT - true iff the IXFAC array is needed
C   --   DON2B - OUT - true iff the IN2ELB array is needed
C   --   DOELED - OUT - true iff the deformed element quarilateral centers
C   --      are needed
C   --   DOELEU - OUT - true iff the undeformed element quarilateral centers
C   --      are needed
C   --   DODEAD - OUT - true iff dead nodes are needed
C   --   DONPS - OUT - true iff node set information is needed
C   --   DOESS - OUT - true iff side set information is needed
C   --   DOSCAL - OUT - true iff the zoom window limits need to be calculated
C   --   MINMSH, MAXMSH - OUT - the minimum and maximum mesh line types
C   --      to be displayed
C   --   MAXHID - OUT - the maximum hidden line option
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DFAC of /DEFORM/
C   --   Uses MSHDEF, MSHNUM, MSHLIN, IHIDOP, NALVAR, DEADNP of /MSHOPT/
C   --   Uses MSCTYP of /MSHLIM/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'mshlim.blk'

      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL

      INTEGER NUMMOD, NDEFVW, IXVW

C   --Set up calculation flags

      ANYDEF = (NUMMOD (MSHDEF, ' ', 'DEFORM', ' ') .GE. 1)
      ANYUND = (NUMMOD (MSHDEF, ' ', 'UNDEFORM', ' ') .GE. 1)
      IF (DFAC .EQ. 0.0) THEN
         ANYUND = .TRUE.
         ANYDEF = .FALSE.
      END IF

      DOIXF = .FALSE.

      DON2B = (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'NODE') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ALL') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'NODE') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ALL') .GE. 1)

      DOELED = (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ELEMENT') .GE. 1)
     &   .OR.  (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ALL') .GE. 1)
      DOELEU = (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ELEMENT') .GE. 1)
     &   .OR.  (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ALL') .GE. 1)
      IF (.NOT. ANYDEF) THEN
         IF (DOELED) DOELEU = .TRUE.
         DOELED = .FALSE.
      END IF

      DODEAD = (NALVAR .GT. 0) .AND. DEADNP

      DONPS = .FALSE.
      DOESS = .FALSE.
      DO 100 IVW = 1, NDEFVW (.FALSE.)
         IVIEW = IXVW (.FALSE., IVW)
         DONPS = DONPS .OR. (NNPSET(IVIEW) .GT. 0)
         DOESS = DOESS .OR. (NESSET(IVIEW) .GT. 0)
  100 CONTINUE

      DOSCAL = (MSCTYP .EQ. 'EACH')

      IF (IS3DIM) THEN
         MAXHID = IHIDOP
         MINMSH = 999
         MAXMSH = 0
         DO 110 IVW = 1, NDEFVW (.FALSE.)
            IVIEW = IXVW (.FALSE., IVW)
            MINMSH = MIN (MINMSH, MSHLIN(IVIEW))
            MAXMSH = MAX (MAXMSH, MSHLIN(IVIEW))
  110    CONTINUE
      END IF

      RETURN
      END
