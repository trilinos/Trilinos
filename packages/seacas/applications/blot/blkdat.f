C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      BLOCK DATA BLKDAT
C=======================================================================

C   --*** BLKDAT *** (BLOT) Block data
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --BLKDAT contains all the COMMON areas defined for BLOT.

C=======================================================================
C                               B L O T
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
C      --These parameters define the indices of 2D and 3D limit arrays

      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'dbnams.blk'
      include 'legopt.blk'
      include 'times.blk'
      include 'layout.blk'
      include 'plcolr.blk'
      include 'cmap-lst.blk'
      include 'light.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

C=======================================================================
C                        M E S H   P L O T

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)
C      --These parameters define the mesh display (see MSHLIN of /MSHOPT/)

      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'cutopt.blk'
      include 'layoud.blk'
      include 'devdat.blk'
      include 'pick.blk'
      include 'sizes.blk'
      include 'linthc.blk'
      include 'sphele.blk'

C=======================================================================
C                             D E T O U R

      include 'detopt.blk'
      include 'cntr.blk'
      include 'icexct.blk'
      include 'icrnbw.blk'
      include 'etcopt.blk'
      include 'nodzom.blk'
      include 'axsplt.blk'
C=======================================================================
C                          X - Y   P L O T

      include 'xyopt.blk'
      include 'xylim.blk'
      include 'xylab.blk'
      include 'neutr.blk'

C=======================================================================
C                        P A T H L I N E

      include 'lnvars.blk'
C=======================================================================
C                             T P L O T

      include 'tpvars.blk'
C=======================================================================
C                             S P L O T

      include 'selne.blk'
      include 'spvars.blk'

C=======================================================================

      DATA IEXCON /0/
      END
