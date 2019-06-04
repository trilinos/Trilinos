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

C $Log: calmag.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1996/06/21 16:06:58  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:55:11  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:49  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CALMAG (DEFOK, IS3DIM, DEFFAC)
C=======================================================================

C   --*** CALMAG *** (MESH) Calculate displacement magnification
C   --   Written by Amy Gilkey - revised 01/26/88
C   --              D. P. Flanagan, 11/17/82
C   --
C   --CALMAG calculates a default magnification factor for the
C   --displacements.
C   --
C   --Parameters:
C   --   DEFOK  - IN  - true iff the displacement variables were found
C   --   IS3DIM - IN  - true if and only if 3D versus 2D mesh
C   --   DEFFAC - OUT - the default displacement magnification
C   --

      LOGICAL DEFOK
      LOGICAL IS3DIM
      REAL DEFFAC

      DEFFAC = 0.0

      IF (.NOT. DEFOK) THEN
         DEFFAC = 0.0
      ELSE
         IF (.NOT. IS3DIM) THEN
            DEFFAC = 1.0
         ELSE
C         --Compute default magnification factor for 3D
            DEFFAC = 1.0
         END IF
      END IF

      RETURN
      END

C*******************************************************************
C      SUBROUTINE CALMAG (A, UNMESH, WHOTIM, NPSURF, DEFOK, IXDEF,
C    &                    IYDEF, IZDEF, DEFFAC)
C   --CALMAG calculates a default magnification factor for the
C   --displacements.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   UNMESH - IN - the limits of the undeformed mesh
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NPSURF - IN - the node numbers of the surface nodes
C   --      or mesh boundary nodes (2D)
C   --   DEFOK - IN - true iff the displacement variables were found
C   --   IXDEF, IYDEF, IZDEF - IN - the indices of the displacement variables
C   --   DEFFAC - OUT - the default displacement magnification
C   --
C   --Common Variables:
C   --   Uses NSTEPS of /DBNUMS/
C   --   Uses IS3DIM, NNPSUR, NUMNPF of /D3NUMS/
C     PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
C     include 'dbnums.blk'
C     include 'd3nums.blk'
C     DIMENSION A(*)
C     REAL UNMESH(KFAR)
C     LOGICAL WHOTIM(*)
C     INTEGER NPSURF(*)
C     LOGICAL DEFOK
C     CHARACTER*8 CONVER
C     IF (.NOT. DEFOK) THEN
C        DEFFAC = 0.0
C     ELSE
C        IF (.NOT. IS3DIM) THEN
C         --Scale displacement variables (SCALER is too inefficient for this)
C            CALL MDRSRV ('CMVAR', KVAR, NUMNPF)
C            CALL MDSTAT (NERR, MEM)
C            IF (NERR .GT. 0) GOTO 110
C            DO 100 ISTEP = 1, NSTEPS
C               IF (WHOTIM(ISTEP)) THEN
C                  CALL GTMVAR (A, IXDEF, -999, ISTEP, NUMNPF, A(KVAR))
C                  CALL MINMAX (NUMNPF, A(KVAR), XMIN, XMAX)
C                  IF (ISTEP .EQ. 1) THEN
C                     XDMIN = XMIN
C                     XDMAX = XMAX
C                  ELSE
C                     XDMIN = MIN (XDMIN, XMIN)
C                     XDMAX = MAX (XDMAX, XMAX)
C                  END IF
C                  CALL GTMVAR (A, IYDEF, -999, ISTEP, NUMNPF, A(KVAR))
C                  CALL MINMAX (NUMNPF, A(KVAR), XMIN, XMAX)
C                  IF (ISTEP .EQ. 1) THEN
C                     YDMIN = XMIN
C                     YDMAX = XMAX
C                  ELSE
C                     YDMIN = MIN (YDMIN, XMIN)
C                     YDMAX = MAX (YDMAX, XMAX)
C                  END IF
C               END IF
C  100       CONTINUE
C            CALL MDDEL ('CMVAR')
C         --Compute default magnification factor for 2D
C           DEFFAC = 1.0
C            TLIM = .05 * MAX (
C     &         UNMESH(KRGT)-UNMESH(KLFT),
C     &         UNMESH(KTOP)-UNMESH(KBOT))
C            DLIM = MAX (ABS(XDMIN), ABS(XDMAX), ABS(YDMIN), ABS(YDMAX))
C            IF (DLIM .NE. 0.0) THEN
C               WRITE (CONVER, 10000) MAX (TLIM / DLIM, 1.0)
C               READ (CONVER, 10000) DEFFAC
C10000           FORMAT (E8.1)
C            ELSE
C               DEFFAC = 0.0
C            END IF
C        ELSE
C         --Compute default magnification factor for 3D
C           DEFFAC = 1.0
C        END IF
C     END IF
C  100 CONTINUE
C     RETURN
C     END
