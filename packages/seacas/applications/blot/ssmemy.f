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

C=======================================================================
      SUBROUTINE SSMEMY (IA,
     &  DONPS, KIDNS, KNNNS, KIXNNS, KLTNNS,
     &  DOESS, KIDSS, KNESS, KNNSS, KIXESS, KIXNSS, KLTESS, KLTNSS)
C=======================================================================

C   --*** SSMEMY *** (MESH) Gets node set and side set info
C   --   Written by Amy Gilkey - revised 04/11/88
C   --   D. P. Flanagan, 11/17/82
C   --
C   --SSMEMY gets the memory for the node sets and side
C   --sets.
C   --
C   --This routine finds the following dynamic memory arrays:
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   DONPS - IN - get node set memory iff true
C   --   KIDNS - OUT - the dynamic memory index of IDNPS
C   --   KNNNS - OUT - the dynamic memory index of NNNPS
C   --   KIXNNS - OUT - the dynamic memory index of IXNNPS
C   --   KLTNNS - OUT - the dynamic memory index of LTNNPS
C   --   DOESS - IN - get side set memory iff true
C   --   KIDSS - OUT - the dynamic memory index of IDESS
C   --   KNESS - OUT - the dynamic memory index of NEESS
C   --   KNNSS - OUT - the dynamic memory index of NNESS
C   --   KIXESS - OUT - the dynamic memory index of IXEESS
C   --   KIXNSS - OUT - the dynamic memory index of IXNESS
C   --   KLTESS - OUT - the dynamic memory index of LTEESS
C   --   KLTNSS - OUT - the dynamic memory index of LTNESS
C   --
C   --Common Variables:
C   --   Uses NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMG/

      include 'dbnumgq.blk'

      DIMENSION IA(*)
      LOGICAL DONPS, DOESS
      logical firstss
      data firstss /.true./

      IF (DONPS) THEN
        CALL MDFIND ('IDNPS', KIDNS, NUMNPS)
        CALL MDFIND ('NNNPS', KNNNS, NUMNPS)
        CALL MDFIND ('IXNNPS', KIXNNS, IDUM)
        CALL MDFIND ('LTNNPS', KLTNNS, N)
      ELSE
        KIDNS = 1
        KNNNS = 1
        KIXNNS = 1
        KLTNNS = 1
      END IF

      IF (DOESS) THEN
        if (firstss) then
          firstss = .false.
          call getssn(ia, ierr)
        end if
        CALL MDFIND ('IDESS', KIDSS, NUMESS)
        CALL MDFIND ('NEESS', KNESS, NUMESS)
        CALL MDFIND ('NNESS', KNNSS, NUMESS)
        CALL MDFIND ('IXEESS', KIXESS, IDUM)
        CALL MDFIND ('IXNESS', KIXNSS, IDUM)
        CALL MDFIND ('LTEESS', KLTESS, N)
        CALL MDFIND ('LTNESS', KLTNSS, N)
      ELSE
        KIDSS = 1
        KNESS = 1
        KNNSS = 1
        KIXESS = 1
        KIXNSS = 1
        KLTESS = 1
        KLTNSS = 1
      END IF

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

 100  CONTINUE
      RETURN
      END
