C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

C=======================================================================
      SUBROUTINE MIRSS (IDFRO, IDBCK,
     &  NESUR, NESFRO, NESBCK, LTEES3, LTSSS3,
     *  COMTOP, NAMELB, NUMELB, IDXELB)
C=======================================================================

C   --*** MIRSS *** (GEN3D) Modifies sideset node order to account
C                           for mirroring about axes
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from WRESS written by Amy Gilkey
C   --
C   --Parameters:
C   --   IDFRO  - IN - ids for front surface side sets; (0) = length
C   --   IDBCK  - IN - ids for back surface side sets; (0) = length
C   --   NESUR  - IN - the number of elements in the surface side set
C   --   NESFRO - IN - the elements in the front surface side set
C   --   NESBCK - IN - the elements in the back surface side set
C   --   LTEES3 - IN - the element faces for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NUMESS of /DBNUMS/
C   --   Uses LESSEO of /DBNUM3/
C   --
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'exodusII.inc'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER NESFRO(*), NESBCK(*)
      INTEGER LTSSS3(*), LTEES3(*)
      CHARACTER*(*)   COMTOP
      CHARACTER*(MXSTLN) NAMELB(NELBLK)
      INTEGER NUMELB(*)
      INTEGER IDXELB(0:*)

      LOGICAL ANYESS

      INTEGER NEWFAC(6)
      INTEGER NEWWED(5)
      DATA NEWFAC /4, 3, 2, 1, 5, 6/
      DATA NEWWED /3, 2, 1, 4, 5/

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYESS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMESS .GT. 0)

C   --Write 3D

C ... If topology is anything other than QUAD or TRI (HEX or WEDGE)
C     it is ignored in this routine.

c ... NOTE: COMTOP is in the 2D topology, NAMELB is in the 3D topology

      IF (ANYESS) THEN
        if (comtop(:8) .eq. 'MULTIPLE') THEN
C ... Sidesets are possibly mixed topology (HEX and WEDGE)
C     Need to iterate the sideset element/side and determine
C     which block the element/side belongs to and then
C     adjust accordingly.
C     This could probably be optimized...

C ... First, set up the index structure to map an element id to the block it is in.
          idxelb(0) = 0
          do i=1, nelblk
            idxelb(i) = numelb(i)
          end do
          do i=1, nelblk
            idxelb(i) = idxelb(i-1) + idxelb(i)
          end do

C ... Now iterate the sideset element/faces
          do nl = 1, lesseo
            iel = ltees3(nl)
            do i=1, nelblk
              if (iel .gt. idxelb(i-1) .and. iel .le. idxelb(i)) then
                iblk = i
                go to 20
              end if
            end do
 20         continue

            if (namelb(iblk)(:1) .eq. 'H') then
              LTSSS3(NL) = NEWFAC(LTSSS3(NL))
            else if (namelb(iblk)(:1) .eq. 'W') then
              LTSSS3(NL) = NEWWED(LTSSS3(NL))
            end if
          end do
        else
C     ... comtop was not equal to "MULTIPLE_TOPOLOGIES", so
C         at this point, the underlying element topology for all elements
C         is the same....

          if (comtop(:4) .eq. 'QUAD') then
C ...       Quad -> hex
            DO NL = 1, LESSEO
C ...         non-front and non-back sidesets
C ...         Front and back don't get mirrored...(?)
              LTSSS3(NL) = NEWFAC(LTSSS3(NL))
            END DO

          else if (comtop(:3) .eq. 'TRI') then
C ...       Tri -> wedge
            DO NL = 1, LESSEO
              LTSSS3(NL) = NEWWED(LTSSS3(NL))
            END DO
          end if
        END IF
      end if
      RETURN
      END
