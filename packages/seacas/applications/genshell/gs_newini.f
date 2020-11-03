C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWINI (IDNSUR, IDESUR, NSSUR, NUMATR)
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
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses NUMNP3, NUMEL3, LNPSN3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses LNPSNO, LESSEO, LESSNO of /DBNUM3/
C   --   Sets NUMNP3, NDIM3, NUMEL3, NELBL3,
C   --      NNPS3, LNPSN3, NESS3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      include 'exodusII.inc'
      INCLUDE 'gs_dbtitl.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      INTEGER NUMATR(NELBLK)

C   --Database title - unchanged

      CONTINUE

C   --Number of dimensions

      NDIM3 = 3

C   --Number of nodes and elements - unchanged

      NUMNP3 = NUMNP
      NUMEL3 = NUMEL

C   --Number of element blocks

      NELBL3 = NELBLK

C   --Lengths of node sets set by NEWNPS
C   --Lengths of side sets set by NEWESS

C   --Number and lengths of sets, including front and back sets

      NNPS3 = NUMNPS + IDNSUR
      LNPSN3 = LNPSNO + IDNSUR*NUMNP
      NESS3 = NUMESS + IDESUR
      LESSE3 = LESSEO + IDESUR*NUMEL
      LESSN3 = LESSNO + IDESUR*NSSUR

C   --Number of attributes per block = 1

      DO 10 IBLK = 1, NELBLK
         NUMATR(IBLK) = 1
 10   CONTINUE

      RETURN
      END
