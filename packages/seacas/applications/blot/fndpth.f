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

C $Log: fndpth.f,v $
C Revision 1.4  2009/03/25 12:36:44  gdsjaar
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
C Revision 1.2  1998/09/03 14:29:52  rwstotz
C fixed problem with the elements path command
C
C Revision 1.1  1994/04/07 20:01:11  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:47  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FNDPTH (NODVAR, LENE, NLNKE, LINKE, XNE, YNE, ZNE,
     &   NPINFO, MAXNE, NNENUM, NENUM)
C=======================================================================

C   --*** FNDPTH *** (SPLOT) Find path between nodes/elements
C   --   Written by Amy Gilkey - revised 10/29/87
C   --
C   --FNDPTH finds a path between nodes or elements.  The nodes
C   --or elements selected are the connected nodes/elements closest to the
C   --straight line between the end points.  "Connected" nodes are in
C   --the same element; "connected" elements share at lease one node.
C   --
C   --Parameters:
C   --   NODVAR - IN - true if nodal versus element plot
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the element connectivity
C   --   XNE, YNE, ZNE - IN - the coordinates of the nodes (if nodes)
C   --      or element centroids (if elements)
C   --   NPINFO - IN/OUT - optimization information
C   --   MAXNE - IN - the number of nodes or elements
C   --   NNENUM - IN/OUT - the number of selected node/element numbers
C   --   NENUM - IN/OUT - the selected node/element numbers
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      LOGICAL NODVAR
      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(NELBLK)
      REAL XNE(*), YNE(*), ZNE(*)
      INTEGER NPINFO(2,*)
      INTEGER NENUM(*)

      INTEGER NECONN(-1:50)
      LOGICAL FIRST
      SAVE FIRST

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
         FIRST = .FALSE.

C      --Calculate the starting and ending element which contains each node

         CALL INIINT (2*NUMNP, 0, NPINFO)
         DO 130 IELB = 1, NELBLK
            IXL0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
            DO 120 IEL = LENE(IELB-1)+1, LENE(IELB)
               IF (LINKE(IXL0+1) .EQ. 0) GOTO 110
               DO 100 K = 1, NLNKE(IELB)
                  N = LINKE(IXL0+K)
                  IF (NPINFO(1,N) .EQ. 0) NPINFO(1,N) = IEL
                  NPINFO(2,N) = IEL
  100          CONTINUE
  110          CONTINUE
               IXL0 = IXL0 + NLNKE(IELB)
  120       CONTINUE
  130    CONTINUE
      END IF

      IF (NNENUM .LT. 2) THEN
         CALL PRTERR ('CMDERR',
     &      'Path must be defined by at least two points')
         GOTO 220
      END IF

C   --Save the path defining points at the end of the list

      IXEND = MAXNE + 1
      DO 140 I = NNENUM, 1, -1
         IXEND = IXEND-1
         NENUM(IXEND) = NENUM(I)
  140 CONTINUE

      ILAST = 0
      NNENUM = 1
      IEND = NENUM(IXEND)
      DO 210 IXEND = IXEND+1, MAXNE
         ISTART = IEND
         IEND = NENUM(IXEND)

C      --Calculate the vector components of start-to-end vector

         UX = XNE(IEND) - XNE(ISTART)
         UY = YNE(IEND) - YNE(ISTART)
         IF (NDIM .EQ. 3) UZ = ZNE(IEND) - ZNE(ISTART)
         IF (NDIM .EQ. 3) THEN
            VECLEN = SQRT (UX*UX + UY*UY + UZ*UZ)
         ELSE
            VECLEN = SQRT (UX*UX + UY*UY)
         END IF
         IF (VECLEN .LE. 0.0) THEN
            CALL PRTERR ('CMDERR',
     &         'Start point is equal to end point')
            GOTO 220
         END IF
         UX = UX / VECLEN
         UY = UY / VECLEN
         IF (NDIM .EQ. 3) UZ = UZ / VECLEN

         SLAST = 0.0
         ITHIS = ISTART

  150    CONTINUE
         IF (ITHIS .NE. IEND) THEN
            IF (NNENUM .GE. IXEND) THEN
               CALL PRTERR ('CMDERR',
     &            'Too many node/element numbers selected')
               GOTO 220
            END IF
            NNENUM = NNENUM + 1

C         --Find all nodes/elements connected to current node/element,
C         --except last node/element

            NCONN = 0
            IF (NODVAR) THEN
               NECONN(0) = ILAST
               DO 170 IEL = NPINFO(1,ITHIS), NPINFO(2,ITHIS)
                  IELB = 0
                  IXL0 = IDBLNK (IELB, IEL, LENE, NLNKE) - 1
                  IL = LOCINT (ITHIS, NLNKE(IELB), LINKE(IXL0+1))
                  IF (IL .GT. 0) THEN
                     DO 160 K = 1, NLNKE(IELB)
                        IF (LINKE(IXL0+K) .EQ. ITHIS) GOTO 160
                        IF (LOCINT (LINKE(IXL0+K), NCONN+1, NECONN(0))
     &                     .LE. 0) THEN
                           NCONN = NCONN + 1
                           NECONN(NCONN) = LINKE(IXL0+K)
                        END IF
  160                CONTINUE
                  END IF
  170          CONTINUE
            ELSE
               NECONN(-1) = ITHIS
               NECONN(0) = ILAST
               ITELB = 0
               IXT0 = IDBLNK (ITELB, ITHIS, LENE, NLNKE) - 1
               DO 190 ILINK = 1, NLNKE(ITELB)
                  INE = LINKE(IXT0+ILINK)
                  DO 180 IEL = NPINFO(1,INE), NPINFO(2,INE)
                     IELB = 0
                     IXL0 = IDBLNK (IELB, IEL, LENE, NLNKE) - 1
                     IL = LOCINT (INE, NLNKE(IELB), LINKE(IXL0+1))
                     IF (IL .GT. 0) THEN
                        IF (LOCINT (IEL, NCONN+2, NECONN(-1)) .LE. 0)
     &                     THEN
                           NCONN = NCONN + 1
                           NECONN(NCONN) = IEL
                        END IF
                     END IF
  180             CONTINUE
  190          CONTINUE
            END IF

C         --Find the node/element with the smallest distance to the
C         --start-to-end vector (DL) that moves toward the end point
C         --(SLAST is distance from start point to where vector from
C         --the last point to start-to-end vector intersect)

            ILAST = ITHIS
            DLMIN = 2 * VECLEN
            DO 200 IX = 1, NCONN
               INE = NECONN(IX)
               XD = XNE(INE) - XNE(ISTART)
               YD = YNE(INE) - YNE(ISTART)
               IF (NDIM .EQ. 3) ZD = ZNE(INE) - ZNE(ISTART)
               IF (NDIM .EQ. 3) THEN
                  S =  UX*XD + UY*YD + UZ*ZD
               ELSE
                  S =  UX*XD + UY*YD
               END IF
               IF (S .GT. SLAST) THEN
                  IF (NDIM .EQ. 3) THEN
                     DL = XD*XD + YD*YD + ZD*ZD - S*S
                  ELSE
                     DL = XD*XD + YD*YD - S*S
                  END IF
                  IF (DL .LT. DLMIN) THEN
                     ITHIS = INE
                     SMIN = S
                     DLMIN = DL
                  END IF
               END IF
  200       CONTINUE

C         --Assign best node/element to list

            IF (ILAST .EQ. ITHIS) THEN
               NENUM(NNENUM) = IEND
               CALL PRTERR ('CMDERR',
     &            'No path found from start to end')
               GOTO 220
            END IF
            SLAST = SMIN
            NENUM(NNENUM) = ITHIS
            GOTO 150
         END IF
  210 CONTINUE

  220 CONTINUE
      RETURN
      END
