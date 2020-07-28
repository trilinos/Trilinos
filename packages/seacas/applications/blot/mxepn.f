C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE  MXEPN (NUME, NLNKE, LINKE, NFPN)
C=======================================================================

C   --*** MXEPN *** (MESH) Count number of elements per node
C   --   Written by Amy Gilkey - revised 10/26/87
C   --
C   --MXEPN counts the maximum number of elements that a node is in.
C   --
C   --Parameters:
C   --   NUME - IN - the number of elements
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity (some nodes may be 0)
C   --   NFPN - SCRATCH - size = 1+NUMNP
C   --
C   --Common Variables:
C   --   Uses NUMNP of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LINKE(NLNKE,NUME)
      INTEGER NFPN(0:NUMNP)

      DO 110 IEL = 1, NUME
         DO 100 ILINK = 1, NLNKE
            INP = LINKE(ILINK,IEL)
            NFPN(INP) = NFPN(INP) + 1
  100    CONTINUE
  110 CONTINUE
      RETURN
      END

      integer function mxepnmx(nfpn)
      include 'dbnums.blk'

      integer nfpn(0:numnp)
      MXEPNMX = 0
      DO 120 INP = 1, NUMNP
         MXEPNMX = MAX (MXEPNMX, NFPN(INP))
  120 CONTINUE
      RETURN
      END
