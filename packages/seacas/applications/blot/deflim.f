C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: deflim.f,v $
C Revision 1.4  2009/03/25 12:36:43  gdsjaar
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
C Revision 1.2  1996/06/21 16:07:07  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:59:22  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:23  gdsjaar
c Added RCS Id and Log to all files
c
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
