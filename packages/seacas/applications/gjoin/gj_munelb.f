C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MUNELB (NELBLK, ISTAT, NUMEL,
     &   IDELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, LINKX, ATRIBX, IXEL, IXELB, NELBX,
     $   ISCR, NAMELB, SCRSTR, LLINK, LATRIB)
C=======================================================================

C   --*** MUNELB *** (GJOIN) Compress and rearrange element blocks
C   --   Written by Amy Gilkey - revised 09/29/87
C   --   Modified by Greg Sjaardema, 07/11/90
C   --      Added element block names
C   --
C   --MUNELB processes the element blocks according to the block status.
C   --Blocks may be combined or deleted.
C   --
C   --Parameters:
C   --   NELBLK - IN/OUT - the number of element blocks
C   --   ISTAT - IN - the status of each block:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NUMEL - IN/OUT - the number of elements
C   --   IDELB - IN/OUT - the element block IDs for each block
C   --   NUMELB - IN/OUT - the number of elements in each block
C   --   NUMLNK - IN/OUT - the number of nodes per element in each block
C   --   NUMATR - IN/OUT - the number of attributes in each block
C   --   LINK - IN/OUT - the connectivity for each block
C   --   ATRIB - IN/OUT - the attributes for each block
C   --   LINKX - SCRATCH - sized to hold the new connectivity
C   --   ATRIBX - SCRATCH - sized to hold the new attributes
C   --   IXEL - OUT - the new element number for each element
C   --   IXELB - SCRATCH - size = NELBLK
C   --   NELBX - SCRATCH - size = NELBLK
C   --   ISCR - SCRATCH - size = NELBLK
C   --   NAMELB - IN/OUT - the names of the element blocks
C   --   SCRSTR - SCRATCH - size = size of NAMELB

      include 'exodusII.inc'

      INTEGER ISTAT(*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*), LINKX(*)
      REAL ATRIB(*), ATRIBX(*)
      INTEGER IXEL(*)
      INTEGER IXELB(*)
      INTEGER NELBX(*)
      INTEGER ISCR(*)
      character*(MXSTLN) namelb(*), scrstr(*)

      DO 100 I = 1, NUMEL
         IXEL(I) = 0
  100 CONTINUE

      JEL0 = 0
      JLNK = 1
      JATR = 1

      JBLK = 0
      DO 140 IELB = 1, NELBLK
         IF (ISTAT(IELB) .EQ. 0) THEN
            NINSET = 1
            ISCR(NINSET) = IELB
         ELSE IF (ISTAT(IELB) .EQ. IELB) THEN
            CALL GETALL (IELB, NELBLK, ISTAT, NINSET, ISCR)
         ELSE
            NINSET = 0
         END IF

         IF (NINSET .GT. 0) THEN
            JBLK = JBLK + 1
            IXELB(JBLK) = IELB
            NELBX(JBLK) = 0
         END IF
         DO 130 ISET = 1, NINSET
            IBLK = ISCR(ISET)
            IEL0 = 0
            ILNK = 1
            IATR = 1
            DO 110 I = 1, IBLK-1
               IEL0 = IEL0 + NUMELB(I)
               ILNK = ILNK + NUMLNK(I) * NUMELB(I)
               IATR = IATR + NUMATR(I) * NUMELB(I)
  110       CONTINUE
            DO 120 I = 1, NUMELB(IBLK)
               IXEL(IEL0+I) = JEL0 + I
  120       CONTINUE
            CALL MOVINT (NUMELB(IBLK) * NUMLNK(IBLK),
     &         LINK(ILNK), LINKX(JLNK))
            CALL MOVREA (NUMELB(IBLK) * NUMATR(IBLK),
     &         ATRIB(IATR), ATRIBX(JATR))
            NELBX(JBLK) = NELBX(JBLK) + NUMELB(IBLK)
            JEL0 = JEL0 + NUMELB(IBLK)
            JLNK = JLNK + NUMLNK(IBLK) * NUMELB(IBLK)
            JATR = JATR + NUMATR(IBLK) * NUMELB(IBLK)
  130    CONTINUE
  140 CONTINUE

      CALL ORDIX  (JBLK, IXELB, NELBLK, IDELB, ISCR, IDELB)
      CALL MOVINT (JBLK, NELBX, NUMELB)
      CALL ORDIX  (JBLK, IXELB, NELBLK, NUMLNK, ISCR, NUMLNK)
      CALL ORDIX  (JBLK, IXELB, NELBLK, NUMATR, ISCR, NUMATR)
      CALL ORDSTR (JBLK, IXELB, NELBLK, NAMELB, SCRSTR, NAMELB)
      NELBLK = JBLK
      NUMEL  = JEL0
      LLINK  = JLNK-1
      LATRIB = JATR-1

      CALL MOVINT (LLINK,  LINKX,  LINK)
      CALL MOVREA (LATRIB, ATRIBX, ATRIB)

      RETURN
      END
