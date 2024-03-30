C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FNDSEL (NODVAR, NNENUM, NENUM,
     &   LENF, NLNKF, LINKF, IF2EL, IX2NP,
     &   NPSURF, NSEL, NPSEL)
C=======================================================================

C   --*** FNDSEL *** (MESH) Find selected nodes/elements
C   --   Written by Amy Gilkey - revised 03/31/88
C   --
C   --FNDSEL constructs a logical array of the selected nodes and elements
C   --that are on the surface of a 3D mesh.
C   --
C   --Parameters:
C   --   NODVAR - IN - true if nodal versus element plot
C   --   NNENUM - IN - the number of selected nodes/elements
C   --   NENUM - IN - the selected nodes/elements
C   --   LENF - IN - the cumulative face counts by element block (if element)
C   --   NLNKF - IN - the number of nodes per face (if element)
C   --   LINKF - IN - the connectivity for all faces (if element)
C   --   IF2EL - IN - the element number of each face (if element)
C   --   IX2NP - IN - the node number for each mesh index (if nodal)
C   --   NPSURF - IN - the node numbers of the surface nodes (if nodal)
C   --   NSEL - OUT - the number of selected nodes
C   --   NPSEL - OUT - the node numbers of the selected nodes
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses NNPSUR, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      include 'd3nums.blk'

      LOGICAL NODVAR
      INTEGER NENUM(NNENUM)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)
      INTEGER IX2NP(NUMNPF)
      INTEGER NPSURF(*)
      INTEGER NPSEL(*)

      INTEGER IFACES(10)

      IF (NODVAR) THEN
         IF (.NOT. IS3DIM) THEN
            CALL INIINT (NUMNPF, 0, NPSEL)
            DO 100 NNP = 1, NNENUM
               INP = NENUM(NNP)
               IF (INP .GT. 0) NPSEL(INP) = 1
  100       CONTINUE

         ELSE
            CALL INIINT (NUMNPF, -1, NPSEL)
            DO 110 IX = 1, NNPSUR
               INP = NPSURF(IX)
               NPSEL(INP) = 0
  110       CONTINUE
            DO 120 NNP = 1, NNENUM
               INP = NENUM(NNP)
               IF (INP .GT. 0) THEN
                  IF (NPSEL(INP) .EQ. 0) NPSEL(INP) = 1
               END IF
  120       CONTINUE
         END IF

      ELSE
         CALL INIINT (NUMNPF, 0, NPSEL)

         DO 150 NEL = 1, NNENUM
            IEL = NENUM(NEL)
            CALL FNDE2F (IEL, LENF, IF2EL, NFARY, IFACES, IELB)
            DO 140 N = 1, NFARY
               IFAC = IFACES(N)
               IELB = 0
               IXL0 = IDBLNK (IELB, IFAC, LENF, NLNKF) - 1
               DO 130 K = 1, NLNKF(IELB)
                  INP = LINKF(IXL0+K)
                  NPSEL(INP) = 1
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
      END IF

C   --Make up the list of selected nodes

      NSEL = 0
      DO 160 INP = 1, NUMNPF
         IF (NPSEL(INP) .GT. 0) THEN
            NSEL = NSEL + 1
            NPSEL(NSEL) = INP
         END IF
  160 CONTINUE

      RETURN
      END
