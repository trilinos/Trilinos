C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MRKNOD (F, HIDENP, XN, YN, ZN, IN2ELB,
     &   NMIN, NMAX, FMIN, FMAX, BLKCOL, *)
C=======================================================================

C   --*** MRKNOD *** (DETOUR) Mark minimum/maximum nodal values
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --MRKNOD marks the nodes which have the minimum or maximum values
C   --with symbols.
C   --
C   --Parameters:
C   --   F - IN - the nodal variable values
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   NMIN, NMAX - IN - the number of variables values matching the
C   --      minimum and the maximum
C   --   FMIN, FMAX - IN - the minimum and maximum function values
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses MAXMIN, MAXMAX of /CNTR/

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'cntr.blk'

      REAL F(NUMNPF)
      LOGICAL HIDENP(NUMNPF)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER IN2ELB(NUMNPF)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT
      LOGICAL MRKMIN, MRKMAX

      CALL UGRCOL (0, BLKCOL)

      MRKMIN = (MAXMIN .GE. NMIN)
      MRKMAX = (MAXMAX .GE. NMAX)

      IF ((FMIN .NE. FMAX) .AND. (MRKMIN .OR. MRKMAX)) THEN

C      --Mark minimum and maximum values with centered symbols

         DO 100 INP = 1, NUMNPF
            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 100
            END IF

            IF (IN2ELB(INP) .GE. 0) THEN
               IF (FMIN .EQ. F(INP)) THEN
                  IF (MRKMIN) THEN
                     IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
                     CALL MPD2SY (1, XN(INP), YN(INP), '\\CCI')
                     CALL MPD2SY (1, XN(INP), YN(INP), '\\CCS')
#else
                     CALL MPD2SY (1, XN(INP), YN(INP), '\CCI')
                     CALL MPD2SY (1, XN(INP), YN(INP), '\CCS')
#endif
                  END IF
               ELSE IF (FMAX .EQ. F(INP)) THEN
                  IF (MRKMAX) THEN
                     IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
                     CALL MPD2SY (1, XN(INP), YN(INP), '\\CX')
                     CALL MPD2SY (1, XN(INP), YN(INP), '\\CCS')
#else
                     CALL MPD2SY (1, XN(INP), YN(INP), '\CX')
                     CALL MPD2SY (1, XN(INP), YN(INP), '\CCS')
#endif
                  END IF
               END IF
            END IF
  100    CONTINUE

         CALL PLTFLU
      END IF

      RETURN
      END
