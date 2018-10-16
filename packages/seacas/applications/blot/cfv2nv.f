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

C $Log: cfv2nv.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
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
C Revision 1.1  1994/04/07 19:55:23  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:57  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CFV2NV (LENF, NLNKF, LINKF, IELBST, ISEVOK, VARFAC,
     &   DOVN2B, IVN2B, IFVCNT, VARNP)
C=======================================================================

C   --*** CFV2NV *** (DETOUR) Convert face variable to nodal variable
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --CFV2NV converts a face variable to a nodal variable by averaging
C   --the values for all "active" elements containing the node.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable of block j exists iff ISEVOK(j)
C   --   VARFAC - IN - the face variable values
C   --   DOVN2B - IN/OUT - set IVN2B and IFVCNT iff true; set false
C   --   IVN2B - IN/OUT - the element block for each node for this variable;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   IFVCNT - IN/OUT - the connected face counts for each node
C   --      for this variable
C   --   VARNP - OUT - the nodal variable values
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses NUMNPF, IS3DIM of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(*)
      REAL VARFAC(*)
      LOGICAL DOVN2B
      INTEGER IVN2B(NUMNPF)
      INTEGER IFVCNT(NUMNPF)
      REAL VARNP(NUMNPF)

      IF (DOVN2B) THEN
         DOVN2B = .FALSE.

         CALL INIINT (NUMNPF, -999, IVN2B)
         CALL INIINT (NUMNPF, 0, IFVCNT)

         DO 120 IELB = 1, NELBLK
            IF ((IELBST(IELB) .GT. 0) .AND. ISEVOK(IELB)) THEN
               IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
               DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  DO 100 K = 1, NLNKF(IELB)
                     INP = LINKF(IXL0+K)
                     IFVCNT(INP) = IFVCNT(INP) + 1
                     IF (IVN2B(INP) .LT. 0) THEN
                        IVN2B(INP) = IELB
                     ELSE IF (IVN2B(INP) .NE. IELB) THEN
                        IVN2B(INP) = 0
                     END IF
  100             CONTINUE
                  IXL0 = IXL0 + NLNKF(IELB)
  110          CONTINUE
            END IF
  120    CONTINUE
      END IF

      CALL INIREA (NUMNPF, 0.0, VARNP)

      DO 150 IELB = 1, NELBLK
         IF ((IELBST(IELB) .GT. 0) .AND. ISEVOK(IELB)) THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
               DO 130 K = 1, NLNKF(IELB)
                  INP = LINKF(IXL0+K)
                  VARNP(INP) = VARNP(INP) + VARFAC(IFAC)
  130          CONTINUE
               IXL0 = IXL0 + NLNKF(IELB)
  140       CONTINUE
         END IF
  150 CONTINUE

      DO 160 INP = 1, NUMNPF
         IF (IFVCNT(INP) .GT. 0)
     &      VARNP(INP) = VARNP(INP) / IFVCNT(INP)
  160 CONTINUE

      RETURN
      END
