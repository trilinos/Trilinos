C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDEB1 (NTXT, IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, natrdm, *)
C=======================================================================

C   --*** RDEB1 *** (TXTEXO) Read database element block misc.
C   --   Written by Amy Gilkey - revised 09/30/87
C   --
C   --RDEB1 reads the element block connectivity and attribute information
C   --from the text file.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   IELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INTEGER LINK(NUMLNK,*)
      REAL ATRIB(natrdm, *)

      CHARACTER*32 STRA, STRB

      NE = 0
      READ (NTXT, *, END=120, ERR=120)
      DO 100 NE = 1, NUMELB
         READ (NTXT, *, END=120, ERR=120) (LINK(I,NE), I=1,NUMLNK)
  100 CONTINUE

      IF (NUMATR .GT. 0) THEN
         NE = 0
         READ (NTXT, *, END=130, ERR=130)
         DO 110 NE = 1, NUMELB
            READ (NTXT, *, END=130, ERR=130) (ATRIB(I,NE), I=1,NUMATR)
  110    CONTINUE
      END IF

      RETURN

  120 CONTINUE
      CALL INTSTR (1, 0, NE, STRA, LSTRA)
      CALL INTSTR (1, 0, IELB, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading NODES for element ' // STRA(:LSTRA)
     &   // ' for BLOCK ' // STRB(:LSTRB))
      GOTO 140
  130 CONTINUE
      CALL INTSTR (1, 0, NE, STRA, LSTRA)
      CALL INTSTR (1, 0, IELB, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading ATTRIBUTES for element ' // STRA(:LSTRA)
     &   // ' for BLOCK ' // STRB(:LSTRB))
      GOTO 140
  140 CONTINUE
      RETURN 1
      END
