C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION GRABRT ()
C=======================================================================

C   --*** GRABRT *** (GRPLIB) Check cancel function (PLT)
C   --   Written by Amy Gilkey - revised 04/28/87
C   --
C   --GRABRT checks the cancel flag.  If it is set, it aborts pending
C   --graphics and returns the terminal to alphanumeric mode.
C   --In any case, the value of the cancel flag is returned as the
C   --function value.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag
C   --   PLTBEL - (PLTLIB) Ring bell
C   --   PLTBGN - (PLTLIB) Erase display surface
C   --   PLTFLU - (PLTLIB) Flush graphics buffer and set alphanumeric mode
C   --   GRSNAP - (GRPLIB) Handle device frame snapping

      LOGICAL CPUIFC

C   --Return cancel flag
      GRABRT = CPUIFC (.FALSE.)

      IF (GRABRT) THEN
C      --Abort plot that is being snapped
         CALL GRSNAP ('ABORT', 0)

C      --Flush buffer (do not erase screen), ring bell, set alphanumeric mode
         CALL PLTFLU
         CALL PLTBEL
         CALL PLTFLU

C      --Send abort message
         WRITE (*, '(1X, A)') '*** Plot set aborted ***'
      END IF

      RETURN
      END
