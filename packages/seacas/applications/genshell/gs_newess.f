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
      SUBROUTINE NEWESS (IDFRO, IDBCK, NSSUR, NLINK, LINK,
     &   ISSFRO, ISSBCK, NSSFRO, NSSBCK)
C=======================================================================

C   $Id: newess.f,v 1.3 1993/11/23 21:19:23 gdsjaar Exp $
C   $Log: newess.f,v $
C   Revision 1.3  1993/11/23 21:19:23  gdsjaar
C   Added capability to handle 8-node quads to 8-node shells
C
c Revision 1.2  1991/03/21  19:37:58  gdsjaar
c Initial change of sidesets from gen3d to genshell
c
c Revision 1.1.1.1  1990/08/20  12:22:16  gdsjaar
c Gen3D Mesh Generation Program
c
c Revision 1.1  90/08/20  12:22:14  gdsjaar
c Initial revision
c

C   --*** NEWESS *** (GEN3D) Calculate 3D side sets
C   --   Written by Amy Gilkey - revised 01/12/88
C   --
C   --NEWESS calculates the side set information for the 3D database.
C   --The elements in the front and back side sets are not calculated
C   --(since they are easily derived), and they are not included in the
C   --length tally.  The nodes in the front and back node sets are
C   --calculated, and their number is returned, but they are not included
C   --in the length tally.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface side sets; (0) = length
C   --   IDBCK - IN - ids for back surface side sets; (0) = length
C   --   LINK - IN - the connectivity for the 2D elements
C   --   ISSFRO - OUT - the elements in the front surface side set
C   --   ISSBCK - OUT - the elements in the back surface side set
C   --   NSSUR - OUT - the number of nodes in the surface side set
C   --   NSSFRO - OUT - the nodes in the front surface side set
C   --   NSSBCK - OUT - the nodes in the back surface side set
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Sets LESSEO, LESSNO of /DBNUM3/
C   --   Uses NNREPL, NEREPL, DIM3 of /PARAMS/

      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER LINK(NLINK,*)
      INTEGER ISSFRO(NUMEL), ISSBCK(NUMEL)
      INTEGER NSSFRO(*), NSSBCK(*)
      INTEGER LTMP(9)

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)

C   --Side set ids - unchanged

      CONTINUE

C   --Side set elements - unchanged

C   --Side set nodes and distribution factors - unchanged

C      --Number of elements in each set

C      --Number of nodes in each set

C   --Number of elements in all sets

C   --Number of nodes in all sets

C   --Set up elements and nodes for front and back side sets

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)

      NSSUR = 0
      IF ((NFRO .GT. 0) .OR. (NBCK .GT. 0)) THEN
        N = 0
        DO 150 JEL = 1, NUMEL
          IF (NFRO .GT. 0) ISSFRO(JEL) = JEL
          IF (NBCK .GT. 0) ISSBCK(JEL) = JEL

          DO 140 I = 1, NLINK
            LTMP(I) = LINK(I,JEL)
 140      CONTINUE
          IF (NFRO .GT. 0) THEN
            NSSFRO(N+1) = LTMP(1)
            NSSFRO(N+2) = LTMP(2)
            NSSFRO(N+3) = LTMP(3)
            NSSFRO(N+4) = LTMP(4)
          END IF
C ... NOTE: Back sidesets must be in reverse order
          IF (NBCK .GT. 0) THEN
            NSSBCK(N+1) = LTMP(4)
            NSSBCK(N+2) = LTMP(3)
            NSSBCK(N+3) = LTMP(2)
            NSSBCK(N+4) = LTMP(1)
          END IF
          IF (NLINK .EQ. 8 .OR. NLINK .EQ. 9) THEN
            IF (NFRO .GT. 0) THEN
              NSSFRO(N+5) = LTMP(5)
              NSSFRO(N+6) = LTMP(6)
              NSSFRO(N+7) = LTMP(7)
              NSSFRO(N+8) = LTMP(8)
            END IF
C ... NOTE: Back sidesets must be in reverse order
            IF (NBCK .GT. 0) THEN
              NSSBCK(N+5) = LTMP(7)
              NSSBCK(N+6) = LTMP(6)
              NSSBCK(N+7) = LTMP(5)
              NSSBCK(N+8) = LTMP(8)
            END IF
            IF (NLINK .EQ. 9) THEN
              IF (NFRO .GT. 0) NSSFRO(N+9) = LTMP(9)
              IF (NBCK .GT. 0) NSSBCK(N+9) = LTMP(9)
            END IF
          END IF
          N = N + NLINK
 150    CONTINUE
        NSSUR = N
      END IF

      RETURN
      END
