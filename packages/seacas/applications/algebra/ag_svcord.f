C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SVCORD (CORD, MAXNE, VARVAL)
C=======================================================================

C   --*** SVCORD *** (ALGEBRA) Save referenced coordinates in VARVAL
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --SVCORD saves the coordinates needed in the equation evaluation
C   --in the VARVAL array.
C   --
C   --Parameters:
C   --   CORD - IN - the coordinates
C   --   MAXNE - IN - the VARVAL dimension (max of NUMEL and NUMNP)
C   --   VARVAL - OUT - the returned coordinates needed
C   --
C   --Common Variables:
C   --   Uses IDVAR, ISTVAR of /VAR../
C   --   Uses NUMNP of /DBNUMS/
C   --   Uses ICOBEG, ICOEND of /DBXVAR/

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'ag_namlen.blk'
      include 'ag_var.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbxvar.blk'

      REAL CORD(NUMNP,NDIM)
      REAL VARVAL(MAXNE,*)

C   --Save any coordinates needed in VARVAL

      DO 100 NVAR = ICOBEG, ICOEND
         ID = IDVAR(NVAR)
         NSTO = ISTVAR(ICURTM,NVAR)
         CALL CPYREA (NUMNP, CORD(1,ID), VARVAL(1,NSTO))
  100 CONTINUE

      RETURN
      END
