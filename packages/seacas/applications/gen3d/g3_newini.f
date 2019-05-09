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
      SUBROUTINE NEWINI (IDNSUR, IDESUR, NSSUR, NESUR, BLKTYP, IBPARM)
C=======================================================================

C   --*** NEWINI *** (GEN3D) Calculate 3D initial variables
C   --   Written by Amy Gilkey - revised 09/02/87
C   --
C   --NEWINI calculates the initial variables for the 3D database.
C   --The output number of nodes and elements and the length of the node
C   --sets and the side sets must be calculated before NEWINI is called.
C   --
C   --Parameters:
C   --   IDNSUR - IN - the number of surface node sets
C   --   IDESUR - IN - the number of surface side sets
C   --   NSSUR - IN - the number of nodes in the surface side set
C   --   NESUR - IN - the number of elements in the surface side set
C   --   BLKTYP - IN - the element block type
C   --   IBPARM - IN - the block parameters (defined by the block type)
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses NUMNP3, NUMEL3, LNPSN3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses LNPSNO, LESSEO, LESSNO of /DBNUM3/
C   --   Sets NUMNP3, NDIM3, NUMEL3, NELBL3,
C   --      NNPS3, LNPSN3, NESS3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      INCLUDE 'g3_dbtitl.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      CHARACTER BLKTYP(NELBLK)
      INTEGER IBPARM(4,NELBLK)

C   --Database title - unchanged

      CONTINUE

C   --Number of dimensions

      NDIM3 = 3

C   --Number of nodes and elements, set by RENUMB

      CONTINUE

C   --Number of element blocks

      NELBL3 = 0
      DO 30 IELB = 1, NELBLK
         IF (BLKTYP(IELB) .EQ. 'T') THEN
            NRSTRT = 1
            NREND = IBPARM(1,IELB) - 1
            NBLK = 1
   10       CONTINUE
            NRSTRT = NREND + 1
            IF (NRSTRT .LE. IBPARM(2,IELB)) THEN
               NREND = MIN (IBPARM(2,IELB), NREND + IBPARM(3,IELB))
            ELSE
               NREND = NEREPL
            END IF
            NBLK = NBLK + 1
            IF (NREND .LT. NEREPL) GOTO 10
            IBPARM(4,IELB) = NBLK
            NELBL3 = NELBL3 + NBLK

         ELSE IF (BLKTYP(IELB) .EQ. 'S') THEN
            NBLK = 1
            NR = NRTRAN(NBLK)
   20       CONTINUE
            IF (NEREPL .GT. NR) THEN
               NBLK = NBLK + 1
               NR = NR + NRTRAN(NBLK)
               GOTO 20
            END IF
            IBPARM(4,IELB) = NBLK
            NELBL3 = NELBL3 + NBLK

         ELSE
            NELBL3 = NELBL3 + 1
         END IF
   30 CONTINUE

C   --Lengths of node sets set by NEWNPS
C   --Lengths of side sets set by NEWESS

C   --LNPSN3 = LNPSNL * NNREPL
C   --LESSE3 = LESSEL * NEREPL
C   --LESSN3 = LESSE3 * 2

C   --Number and lengths of sets, including front and back sets

      NNPS3 = NUMNPS + IDNSUR
      LNPSN3 = LNPSNO + IDNSUR*NUMNP

      NESS3 = NUMESS + IDESUR
      LESSE3 = LESSEO + IDESUR*NESUR
      LESSN3 = LESSNO + IDESUR*NSSUR

      RETURN
      END
