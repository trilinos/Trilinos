C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Log: dtread.f,v $
C Revision 1.3  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 19:59:57  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:39  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE DTREAD (A, ISTEP, IDTVAR, NNDVAR, NEDVAR,
     &   LENF, IF2EL, IELBST, ISEVOK, VARNP, VARFAC, VAR, LVARF)
C=======================================================================

C   --*** DTREAD *** (DETOUR) Read variables for time step
C   --   Written by Amy Gilkey - revised 05/11/88
C   --
C   --DTREAD reads the nodal and element variables needed for the time step.
C   --The element variables are converted to face variables.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the selected time step number
C   --   IDTVAR - IN - the variable numbers
C   --   NNDVAR - IN - the number of nodal variables needed
C   --   NEDVAR - IN - the number of element variables needed
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VARNP - OUT - the nodal variable values
C   --      (sized for all needed nodal variables)
C   --   VARFAC - OUT - the face variable values
C   --      (sized for all needed element variables)
C   --   VAR - SCRATCH - array to hold an element variable
C   --
C   --Common Variables:
C   --   Uses NUMNP, NUMEL, NELBLK of /DBNUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      DIMENSION A(*)
      INTEGER IDTVAR(*)
      INTEGER LENF(0:NELBLK)
      INTEGER IF2EL(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      REAL VARNP(NUMNP,*)
      REAL VARFAC(LVARF,*)
      REAL VAR(NUMEL)

      CHARACTER TYP
      INTEGER IX(3)

      DO 160 IVAR = MAX (NNDVAR, NEDVAR), 1, -1

C      --Get the variable type

         CALL DBVTYP_BL (IDTVAR(IVAR), TYP, ID)

C      --Read in variable

         IF (TYP .EQ. ' ') THEN

C         --Zero out nodal variable and element variable if needed

            IF (NNDVAR .GE. IVAR) THEN
               DO 110 INP = 1, NUMNP
                  VARNP(INP,IVAR) = 0.0
  110          CONTINUE
            END IF

            IF (NEDVAR .GE. IVAR) THEN
               DO 130 IELB = 1, NELBLK
                  IF (IELBST(IELB) .GT. 0) THEN
                     DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
                        VARFAC(IFAC,IVAR) = 0.0
  120                CONTINUE
                  END IF
  130          CONTINUE
            END IF

         ELSE IF (TYP .EQ. 'N') THEN

C         --Get nodal variable

            CALL GTMVAR (A, IDTVAR(IVAR), -999, ISTEP, NUMNPF,
     &         VARNP(1,IVAR))

         ELSE IF (TYP .EQ. 'E') THEN

C            --Get element variable, change to face variable

               CALL GETVAR (A, IDTVAR(IVAR), -1, ISTEP, NUMEL, VAR)

               DO 150 IELB = 1, NELBLK
                  IF ((IELBST(IELB) .GT. 0) .AND. ISEVOK(IELB,ID)) THEN
                     DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
                        VARFAC(IFAC,IVAR) = VAR(IF2EL(IFAC))
  140                CONTINUE
                  END IF
  150          CONTINUE

         END IF
  160 CONTINUE

      RETURN
      END
