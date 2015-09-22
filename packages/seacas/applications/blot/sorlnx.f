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

C $Log: sorlnx.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:13:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:07  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SORLNX (NSETS, NEWELB, IE2ELB,
     &   NLNKF, LENSC, LINKSC, IF2ESC, IF22SC,
     &   LENF, LINKF, IF2EL, IF2EL2)
C=======================================================================

C   --*** SORLNX *** (MESH) Sort faces
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --SORLNX resorts all the faces into new element block sets.
C   --
C   --Parameters:
C   --   NSETS - IN - the number of element block sets in LINKF
C   --   NEWELB - IN - the new element block set for each face
C   --   IE2ELB - IN - the element for each element block
C   --   NLNKF - IN - the number of nodes for each element block
C   --   LENSC - IN - the cumulative face counts by element block for LINKSC
C   --   LINKSC - IN - the unsorted connectivity for all faces
C   --   IF2ESC - IN - the element number of each face in LINKSC
C   --   IF22SC - IN - the secondary element number of each face in LINKSC
C   --   LENF - OUT - the cumulative face counts by element block
C   --   LINKF - OUT - the connectivity for all faces
C   --   IF2EL - OUT - the element number of each face
C   --   IF2EL2 - OUT - the secondary element number of each face
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER NEWELB(*)
      INTEGER IE2ELB(*)
      INTEGER NLNKF(NELBLK)
      INTEGER LENSC(0:NSETS)
      INTEGER LINKSC(*)
      INTEGER IF2ESC(*), IF22SC(*)
      INTEGER LENF(0:NSETS)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)

C   --Move faces into appropriate element block

      IX = 0
      IXL0 = 0
      DO 120 NELB = 1, NSETS

         IXLSC0 = 0
         DO 110 IELB = 1, NSETS
            IF (IELB .LE. NELBLK) NL = NLNKF(IELB)

            DO 100 IFAC = LENSC(IELB-1)+1, LENSC(IELB)
               IF (IELB .GT. NELBLK) NL = NLNKF(IE2ELB(IF2ESC(IFAC)))

               IF (NELB .EQ. IABS (NEWELB(IFAC))) THEN

                  IF (NEWELB(IFAC) .LT. 0) THEN

C                  --IF2ESC(IFAC) holds the "live" element number
                     I = IF2ESC(IFAC)
                     IF2ESC(IFAC) = IF22SC(IFAC)
                     IF22SC(IFAC) = I

C                  --Swap nodes to simulate surface being defined
C                  --by facing element
                     IF (NL .EQ. 4) THEN
                        I = LINKSC(IXLSC0+2)
                        LINKSC(IXLSC0+2) = LINKSC(IXLSC0+4)
                        LINKSC(IXLSC0+4) = I
                     END IF
                  END IF

                  IX = IX + 1
                  CALL CPYINT (NL, LINKSC(IXLSC0+1), LINKF(IXL0+1))
                  IXL0 = IXL0 + NL
                  IF2EL(IX) = IF2ESC(IFAC)
                  IF (IS3DIM) IF2EL2(IX) = IF22SC(IFAC)
               END IF

               IXLSC0 = IXLSC0 + NL
  100       CONTINUE

  110    CONTINUE
         LENF(NELB) = IX
  120 CONTINUE
      RETURN
      END
