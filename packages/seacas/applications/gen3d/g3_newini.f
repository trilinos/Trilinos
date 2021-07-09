C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWINI (IDNSUR, IDESUR, NSSUR, NESUR, BLKTYP, IBPARM)
C=======================================================================

C   --*** NEWINI *** (GEN3D) Calculate 3D initial variables
C   --   Written by Amy Gilkey - revised 09/02/87
C   --
C   --NEWINI calculates the initial variables for the 3D database.
C   --The output number of nodes and elements and the length of the node
C   --sets and the side sets must be calculated before NEWINI is called.
C   --
C   --Parameters:
C   --   IDNSUR - IN - the number of surface node sets
C   --   IDESUR - IN - the number of surface side sets
C   --   NSSUR - IN - the number of nodes in the surface side set
C   --   NESUR - IN - the number of elements in the surface side set
C   --   BLKTYP - IN - the element block type
C   --   IBPARM - IN - the block parameters (defined by the block type)
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses NUMNP3, NUMEL3, LNPSN3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses LNPSNO, LESSEO, LESSNO of /DBNUM3/
C   --   Sets NUMNP3, NDIM3, NUMEL3, NELBL3,
C   --      NNPS3, LNPSN3, NESS3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      INCLUDE 'g3_dbtitl.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      CHARACTER BLKTYP(NELBLK)
      INTEGER IBPARM(4,NELBLK)

C   --Database title - unchanged

      CONTINUE

C   --Number of dimensions

      NDIM3 = 3

C   --Number of nodes and elements, set by RENUMB

      CONTINUE

C   --Number of element blocks

      NELBL3 = 0
      DO 30 IELB = 1, NELBLK
         IF (BLKTYP(IELB) .EQ. 'T') THEN
            NRSTRT = 1
            NREND = IBPARM(1,IELB) - 1
            NBLK = 1
   10       CONTINUE
            NRSTRT = NREND + 1
            IF (NRSTRT .LE. IBPARM(2,IELB)) THEN
               NREND = MIN (IBPARM(2,IELB), NREND + IBPARM(3,IELB))
            ELSE
               NREND = NEREPL
            END IF
            NBLK = NBLK + 1
            IF (NREND .LT. NEREPL) GOTO 10
            IBPARM(4,IELB) = NBLK
            NELBL3 = NELBL3 + NBLK

         ELSE IF (BLKTYP(IELB) .EQ. 'S') THEN
            NBLK = 1
            NR = NRTRAN(NBLK)
   20       CONTINUE
            IF (NEREPL .GT. NR) THEN
               NBLK = NBLK + 1
               NR = NR + NRTRAN(NBLK)
               GOTO 20
            END IF
            IBPARM(4,IELB) = NBLK
            NELBL3 = NELBL3 + NBLK

         ELSE
            NELBL3 = NELBL3 + 1
         END IF
   30 CONTINUE

C   --Lengths of node sets set by NEWNPS
C   --Lengths of side sets set by NEWESS

C   --LNPSN3 = LNPSNL * NNREPL
C   --LESSE3 = LESSEL * NEREPL
C   --LESSN3 = LESSE3 * 2

C   --Number and lengths of sets, including front and back sets

      NNPS3 = NUMNPS + IDNSUR
      LNPSN3 = LNPSNO + IDNSUR*NUMNP

      NESS3 = NUMESS + IDESUR
      LESSE3 = LESSEO + IDESUR*NESUR
      LESSN3 = LESSNO + IDESUR*NSSUR

      RETURN
      END
