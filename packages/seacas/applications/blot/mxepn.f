C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: mxepn.f,v $
C Revision 1.4  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.2  1999/03/09 22:02:50  gdsjaar
C Fixed problem with incorrect 'mesh contiguity' errors. Max number of
C faces connected to a node was being based on block-by-block counting
C of elements connected to a node which didn't work well for meshes
C containing lots of blocks. Changed to do count based on all element blocks.
C
C Revision 1.1  1994/04/07 20:06:04  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:12  gdsjaar
c Added RCS Id and Log to all files
c
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
