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

c
C=======================================================================
      SUBROUTINE SYMBOL_BL (ISYTYP, VARFAC, LENF, NLNKF, HIDEF,
     &   XF, YF, ZF, ISVOK, BLKCOL, IDELB, *)
C=======================================================================

C   --*** SYMBOL *** (DETOUR) Plot element symbol
C   --   Written by Amy Gilkey, revised 10/12/87
C   --   D. P. Flanagan, 12/08/83
C   --
C   --SYMBOL drives the user symbol interface for element variables.
C   --It processes each element by element block, computing scaling factors
C   --and element information, then calling the user symbol routine.
C   --Only elements of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   ISYTYP - IN - the symbol type (as in MODTYP of /DETOPT/)
C   --   VARFAC - IN - the face variable values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   ISVOK - IN - ISVOK(i) is true iff the symbol variables are defined
C   --      for element block i
C   --   BLKCOL - IN - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses VECSCL of /ETCOPT/
C   --   Uses DTW, VWSCL of /DEVDAT/

      PARAMETER (KSCHSZ=2)

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /ETCOPT/ VECSCL
      COMMON /DEVDAT/ DTW, VWSCL

      CHARACTER*(*) ISYTYP
      REAL VARFAC(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      LOGICAL HIDEF(*)
      REAL XF(*), YF(*), ZF(*)
      LOGICAL ISVOK(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      LOGICAL GRABRT
      LOGICAL PLTSTT, PLTGTT, LDUM

      SSCL = VECSCL * DTW

C   --Set symbol size, must be reset before exit
      LDUM = PLTGTT (KSCHSZ, SYMSZ)
      LDUM = PLTSTT (KSCHSZ, VECSCL*SYMSZ)

      DO 110 IELB = 1, NELBLK
         IF (ISVOK(IELB)) THEN

C         --Set the symbol color

c            CALL UGRCOL (IDELB(IELB), BLKCOL)
            CALL UGRCOL (IELB, BLKCOL)

            DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (IS3DIM) THEN
                  IF (HIDEF(IFAC)) GOTO 100
               END IF

               IF (GRABRT ()) GOTO 120

C            --Call symbol routine

               CALL USRSYM (ISYTYP, IS3DIM,
     &            XF(IFAC), YF(IFAC), ZF(IFAC), VARFAC(IFAC), SSCL)

  100       CONTINUE

            CALL PLTFLU
         END IF
  110 CONTINUE

C   --Reset symbol size
      LDUM = PLTSTT (KSCHSZ, SYMSZ)

      RETURN

  120 CONTINUE
C   --Reset symbol size
      LDUM = PLTSTT (KSCHSZ, SYMSZ)

      RETURN 1
      END
