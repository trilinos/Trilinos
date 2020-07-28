C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DEFLIM (A, WHOTIM, XN, YN, ZN, NPSURF)
C=======================================================================

C   --*** DEFLIM *** (MESH) Calculate deformed mesh limits
C   --   Written by Amy Gilkey - revised 01/22/88
C   --   D. P. Flanagan, 11/17/82
C   --
C   --DEFLIM calculates the limits of the deformed mesh.
C   --
C   --Parameters:
C   --   A        - IN - the dynamic memory base array
C   --   WHOTIM   - IN - true iff whole (versus history) time step
C   --   XN,YN,ZN - IN - the nodal coordinates (ZN for 3D only)
C   --   NPSURF   - IN - the node numbers of the surface nodes
C   --                   or mesh boundary nodes (2D)
C   --
C   --Common Variables:
C   --   Uses NDIM, NVARNP, NSTEPS of /DBNUMS/
C   --   Uses IS3DIM, NNPSUR, NUMNPF of /D3NUMS/
C   --   Sets DEFOK, DEFFAC, IXDEF, IYDEF, IZDEF of /DEFORM/
C   --   Sets ALMESH of /MSHLIM/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshlim.blk'

      DIMENSION A(*)
      LOGICAL WHOTIM(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER NPSURF(*)

C   --Calculate the default magnification

      CALL CALMAG (DEFOK, IS3DIM, DEFFAC)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 110

C   --Find deformed mesh limits

C   --Note that the displacement variables will be put on the random file
C   --which will improve efficiency

      CALL CPYREA (2*NDIM, UNMESH, ALMESH)

      IF (DEFOK) THEN
         CALL MDRSRV ('DLX', KDXN, NUMNPF)
         CALL MDRSRV ('DLY', KDYN, NUMNPF)
         IF (IS3DIM) THEN
            CALL MDRSRV ('DLZ', KDZN, NUMNPF)
         ELSE
            KDZN = 1
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 110

         DO 100 ISTEP = 1, NSTEPS
            IF (WHOTIM(ISTEP)) THEN
               CALL DEFXYZ (A, ISTEP, DEFFAC, .TRUE., NPSURF,
     &              XN, YN, ZN, A(KDXN), A(KDYN), A(KDZN))
               CALL MINMXS (NNPSUR, NPSURF, A(KDXN), XMIN, XMAX)
               ALMESH(KLFT) = MIN (ALMESH(KLFT), XMIN)
               ALMESH(KRGT) = MAX (ALMESH(KRGT), XMAX)
               CALL MINMXS (NNPSUR, NPSURF, A(KDYN), XMIN, XMAX)
               ALMESH(KBOT) = MIN (ALMESH(KBOT), XMIN)
               ALMESH(KTOP) = MAX (ALMESH(KTOP), XMAX)
               IF (IS3DIM) THEN
                  CALL MINMXS (NNPSUR, NPSURF, A(KDZN), XMIN, XMAX)
                  ALMESH(KNEA) = MIN (ALMESH(KNEA), XMIN)
                  ALMESH(KFAR) = MAX (ALMESH(KFAR), XMAX)
               END IF
            END IF
  100    CONTINUE

         CALL MDDEL ('DLX')
         CALL MDDEL ('DLY')
         IF (IS3DIM) THEN
            CALL MDDEL ('DLZ')
         END IF
      END IF

  110 CONTINUE
      RETURN
      END
