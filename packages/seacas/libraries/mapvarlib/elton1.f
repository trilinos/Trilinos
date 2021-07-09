C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*     DECK, ELTON1
      SUBROUTINE ELTON1(CNTRA,SOLEA,SOLENA,IDBLK,
     &     XA,YA,ZA,ICONA,NDLSTA,
     &     INVLN,INVCN,MAXLN,ISTP,
     &     ITT, iblk)

C     *********************************************************************

C     Subroutine ELTON1 extracts nodal values of element variables by
C     c  performing a weighted least squares fit (4 or more elements) or
C     c  a triangulation (3 elements) over the centroids of the elements
C     c  attached to the current node.

C     Each element block must be processed independently in order to
C     avoid averaging element variables across material boundaries.
C     Note: the last set of DO loops acts over all nodes; to make sense
C     one element block must be completely processed before another
C     element block is sent into this subroutine.

C     Calls subroutines CNTR, VOL, EXTQ, AVG, EXTH, ERROR

C     Called by MAPVAR

C     *********************************************************************

C     CNTRA     a list of element centroid coordinates for all elements
C     in the current element block (1:ndima,1:numeba)
C     SOLEA     element variables (1:numeba,1:nvarel)
C     SOLENA    element variables at nodes (1:nodesa,1:nvarel)
C     IDBLK     current element block I.D.
C     XA,YA,ZA  coordinates
C     ICONA      mesh-A connectivity (1:nelnda,1:numeba)
C     NDLSTA    list of nodes in element block - from RDA2 (1:numnda)
C     INVLN     number of elements per node (1:numnda)
C     INVCN     inverse connectivity (1:nelnda,1:numnda)
C     MAXLN     maximum number of elements connected to any node
C     ITT       truth table
C     iblk      element block being processed (not ID)

C**   RELATIONSHIP BETWEEN NODAL IDENTIFICATIONS **
C     IGLND  = NDLSTA(INOD)   = ICONA(NOWLND,INVCN(1,IGLND))

C     *********************************************************************

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'

      DIMENSION CNTRA(NUMEBA,*), SOLEA(NUMEBA,*)
      DIMENSION SOLENA(NODESA,NVAREL)
      DIMENSION XX(27), YY(27), ZZ(27), IFCLND(4), IEGLND(2)
      DIMENSION XA(*), YA(*), ZA(*), ICONA(NELNDA,*), NDLSTA(*)
      DIMENSION INVCN(MAXLN,*), INVLN(*), ITT(NVAREL,*)

C     *********************************************************************

      NXTLND = 0
      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
         CALL ERROR('ELTON1','ELEMENT TYPE',' ',ITYPE,
     &        'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &        0,' ',' ',1)
      END IF

      DO I = 1, NODESA
         DO J = 1, NVAREL
            SOLENA(I,J) = 0.
         end do
      end do

C     load up CNTRA array - coordinates of mesh-A element centroids

C     NNODES = NNELM(ITYPE)
         NNODES = NELNDA
         IF (ITYPE .EQ. 6) NNODES = 4
         IF (NDIMA .EQ. 2) THEN
            DO 40 IEL = 1, NUMEBA
               DO 50 I = 1, NNODES
                  INODE = ICONA(I,IEL)
                  XX(I) = XA(INODE)
                  YY(I) = YA(INODE)
                  ZZ(I) = 0.
 50            CONTINUE
               CALL CNTR(ITYPE,XX,YY,ZZ,CNTRA(IEL,1),CNTRA(IEL,2),DUMMY)
 40         CONTINUE
         ELSE
            DO 60 IEL = 1, NUMEBA
               DO 70 I = 1, NNODES
                  INODE = ICONA(I,IEL)
                  XX(I) = XA(INODE)
                  YY(I) = YA(INODE)
                  ZZ(I) = ZA(INODE)
 70            CONTINUE
               CALL CNTR(ITYPE,XX,YY,ZZ,CNTRA(IEL,1),CNTRA(IEL,2),
     &              CNTRA(IEL,3))
 60         CONTINUE
         END IF

C     put element variables into SOLEA array

         DO 80 IVAR = 1, NVAREL
            IF (ITT(IVAR,iblk) .EQ. 0)GO TO 80
            CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)

            IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN

C     replace element mass with density

               DO 90 IEL = 1, NUMEBA
                  DO 100 I = 1, NNODES
                     XX(I) = XA(ICONA(I,IEL))
                     YY(I) = YA(ICONA(I,IEL))
                     IF (NDIMA .EQ. 3)THEN
                        ZZ(I) = ZA(ICONA(I,IEL))
                     ELSE
                        ZZ(I) = 0.
                     END IF
 100              CONTINUE
                  CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
                  SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
 90            CONTINUE
            END IF
 80      CONTINUE

C     start least squares extrapolation

c     First check element type
c     3 = 4-node quad (2d)
c     10 = 8-node hex  (3d)

c*******

         IF (ITYPE .EQ. 3)THEN

C     Find the elements connected to the node. If fewer than 3 elements,
C     adjust search to find additional elements. If unable to get at
C     least 3 elements, must be treated as special case (just average
C     element values at node)(see below).

            DO 110 INOD = 1, NUMNDA
               IGLND = NDLSTA(INOD)

C     Process special case of only 1 element attached to node

               IF (INVLN(IGLND) .EQ. 1)THEN

C     Get node number diagonally across element, in most cases this
C     node will have 4 elements attached.

                  NXTLND = 0
                  DO 120 I = 1, NNODES
                     IF (IGLND .EQ. ICONA(I,INVCN(1,IGLND))) THEN
                        NXTLND = I + 2
                     END IF
 120              CONTINUE
                  IF (NXTLND .GT. NNODES) NXTLND = NXTLND - NNODES
                  NXGLND = ICONA(NXTLND,INVCN(1,IGLND))

C     If 3 or more elements perform least
c     squares extrapolation to original node. If 2 or less elements,
c     average original element variables at original node

                  IF (INVLN(NXGLND) .GT. 2)THEN
                     CALL EXTQ(IGLND,INVCN,MAXLN,NXGLND,INVLN(NXGLND),
     $                    XA,YA, CNTRA,SOLEA,SOLENA,ITT,iblk)
                  ELSE
                     CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                    SOLEA,SOLENA,ITT,iblk)
                  END IF

C     Process special case of only 2 elements attached to node

               ELSE IF (INVLN(IGLND) .EQ. 2)THEN

c     get second node that is shared by both elements. That is the
c     node on the other end of the shared element side.

                  DO I = 1, NNODES
                     DO J = 1, NNODES
                        IF(ICONA(I,INVCN(1,IGLND)) .NE. IGLND .AND.
     &                       ICONA(I,INVCN(1,IGLND)) .EQ.
     $                       ICONA(J,INVCN(2,IGLND))) THEN
                           NXGLND = ICONA(I,INVCN(1,IGLND))
                        END IF
                     end do
                  end do

c     If this second node has more than 2 elements, extrapolate. Otherwise
c     average. (at original node)

                     IF (INVLN(NXGLND) .GT. 2)THEN
                        CALL EXTQ(IGLND,INVCN,MAXLN,NXGLND,
     $                       INVLN(NXGLND), XA,YA,CNTRA,
     $                       SOLEA,SOLENA,ITT,iblk)
                     ELSE
                        CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                       SOLEA,SOLENA,ITT,iblk)
                     END IF
                  ELSE
                     CALL EXTQ(IGLND,INVCN,MAXLN,IGLND,INVLN(IGLND),
     $                    XA,YA,CNTRA,SOLEA,SOLENA,ITT,iblk)
                  END IF
 110           CONTINUE

c*****

            ELSE IF (ITYPE .EQ. 10)THEN

c     Do for 8-node hex in 3D, similar to 4-node quad in 2D above

               DO 200 INOD = 1, NUMNDA
                  IGLND = NDLSTA(INOD)

c     First find elements connected to node - inverse connectivity

C     Similar to 2D, process special cases

                  NOWLND = 0
                  DO 210 I = 1, NNODES
                     IF (IGLND .EQ. ICONA(I,INVCN(1,IGLND)))THEN
                        NOWLND = I
                        GO TO 220
                     END IF
 210              CONTINUE
 220              CONTINUE

C     Only 1 element connected to node, find node diagonally across hex

                  IF (INVLN(IGLND) .EQ. 1)THEN
                     IF (NOWLND .EQ. 1 .OR. NOWLND .EQ. 2)THEN
                        NXTLND = NOWLND + 6
                     ELSE IF (NOWLND .EQ. 3 .OR. NOWLND .EQ. 4)THEN
                        NXTLND = NOWLND + 2
                     ELSE IF (NOWLND .EQ. 5 .OR. NOWLND .EQ. 6)THEN
                        NXTLND = NOWLND - 2
                     ELSE IF (NOWLND .EQ. 7 .OR. NOWLND .EQ. 8)THEN
                        NXTLND = NOWLND - 6
                     END IF
                     NXGLND = ICONA(NXTLND,INVCN(1,IGLND))

                     IF (INVLN(NXGLND) .GT. 5)THEN
                        CALL EXTH(IGLND,INVCN,MAXLN,NXGLND,
     $                       INVLN(NXGLND),XA,YA,ZA,CNTRA,
     $                       SOLEA,SOLENA,ITT,iblk)
                     ELSE
                        CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                       SOLEA,SOLENA,ITT,iblk)
                     END IF
                     go to 200

C     Only 2 elements connected to node, find node diagonally across
C     shared face of 2 elements

                  ELSE IF (INVLN(IGLND) .EQ. 2)THEN
                     DO 250 J = 1, NNODES
                        IF (NOWLND .EQ. 1)THEN
                           IF (ICONA(3,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 3
                              GO TO 260
                           ELSE IF (ICONA(6,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 6
                              GO TO 260
                           ELSE IF (ICONA(8,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 8
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 2)THEN
                           IF (ICONA(4,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 4
                              GO TO 260
                           ELSE IF (ICONA(5,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 5
                              GO TO 260
                           ELSE IF (ICONA(7,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 7
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 3)THEN
                           IF (ICONA(1,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 1
                              GO TO 260
                           ELSE IF (ICONA(6,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 6
                              GO TO 260
                           ELSE IF (ICONA(8,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 8
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 4)THEN
                           IF (ICONA(2,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 2
                              GO TO 260
                           ELSE IF (ICONA(5,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 5
                              GO TO 260
                           ELSE IF (ICONA(7,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 7
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 5)THEN
                           IF (ICONA(2,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 2
                              GO TO 260
                           ELSE IF (ICONA(4,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 4
                              GO TO 260
                           ELSE IF (ICONA(7,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 7
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 6)THEN
                           IF (ICONA(1,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 1
                              GO TO 260
                           ELSE IF (ICONA(3,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 3
                              GO TO 260
                           ELSE IF (ICONA(8,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 8
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 7)THEN
                           IF (ICONA(2,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 2
                              GO TO 260
                           ELSE IF (ICONA(4,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 4
                              GO TO 260
                           ELSE IF (ICONA(5,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 5
                              GO TO 260
                           END IF
                        ELSE IF (NOWLND .EQ. 8)THEN
                           IF (ICONA(1,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 1
                              GO TO 260
                           ELSE IF (ICONA(3,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 3
                              GO TO 260
                           ELSE IF (ICONA(6,INVCN(1,IGLND)) .EQ.
     &                             ICONA(J,INVCN(2,IGLND)))THEN
                              NXTLND = 6
                              GO TO 260
                           END IF
                        ELSE
                           CALL ERROR('ELTON1',
     $                          'IN 3-D, LEAST-SQUARES, 2 ELEMENTS',
     &               'NEXT NODE - DIAGONALLY ACROSS FACE NOT FOUND',
     $                          0,' ',0,' ',' ',1)
                        END IF
 250                 CONTINUE
 260                 CONTINUE
                     NXGLND = ICONA(NXTLND,INVCN(1,IGLND))

                     IF (INVLN(NXGLND) .GT. 5)THEN
                        CALL EXTH(IGLND,INVCN,MAXLN,NXGLND,
     $                       INVLN(NXGLND),XA,YA,ZA,CNTRA,
     $                       SOLEA,SOLENA,ITT,iblk)
                     ELSE
                        CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                       SOLEA,SOLENA,ITT,iblk)
                     END IF
                     go to 200
                  ELSE IF (INVLN(IGLND) .LT. 8)THEN

C     If 3 to 7 elements are connected to a node, check for shared edge.
C     If all elements share an edge, transfer to other end of that edge.
C     Otherwise, extrapolate/average with what you have.

C     Step-1 find shared face or edge between element-1 and element-2

                     K = 0
                     DO 300 I = 1, NNODES
                        DO 310 J = 1, NNODES
                           IF (ICONA(I,INVCN(1,IGLND)) .EQ.
     &                          ICONA(J,INVCN(2,IGLND)))THEN
                              K = K + 1
                              IFCLND(K) = I
                              GO TO 300
                           END IF
 310                    CONTINUE
 300                 CONTINUE

C     If K=4, shared face, process element 3 to determine shared edge
C     If K=6, or K=8, then there is a degenerate hex in the mesh -
C     don't worry how it is processed
C     If K=2, shared edge
C     IF K=1, no shared edge, extrapolate/average with elements you got

                     KC = 0
                     IF (K .EQ. 6 .OR. K .EQ. 8)THEN
                        GO TO 380
                     ELSE IF (K .EQ. 4)THEN
                        DO 320 IC = 1, K
                           DO 330 JC = 1, NNODES
                              IF (ICONA(IFCLND(IC),INVCN(1,IGLND)) .EQ.
     &                             ICONA(JC,INVCN(3,IGLND)))THEN
                                 KC = KC + 1
                                 IEGLND(KC) = IFCLND(IC)
                                 GO TO 320
                              END IF
 330                       CONTINUE
 320                    CONTINUE
                     ELSE IF (K .EQ. 2) THEN
                        KC = 2
                        IEGLND(1) = IFCLND(1)
                        IEGLND(2) = IFCLND(2)
                     ELSE IF (K .EQ. 1)THEN
                        GO TO 380
                     ELSE
                        CALL ERROR('ELTON1',
     $                       'IN 3-D, LEAST-SQUARES, <8 ELEMENTS',
     $             'EDGE CHECKING - K MUST BE EVEN NO., OR 1 K= ',
     $                       K,' ',0,' ',' ',1)
                     END IF
                     DO 340 IL = 1, INVLN(IGLND)
                        IEDGCT = 0
                        DO 350 JJ = 1, NNODES
                          IF ( KC .GE. 2) THEN
                            IF (ICONA(IEGLND(1),INVCN(1,IGLND)) .EQ.
     &                        ICONA(JJ,INVCN(IL,IGLND)) .OR.
     &                        ICONA(IEGLND(2),INVCN(1,IGLND)) .EQ.
     &                        ICONA(JJ,INVCN(IL,IGLND)))THEN
                              IEDGCT = IEDGCT + 1
                            END IF
                          ELSE IF ( KC .EQ. 1) THEN
                            IF (ICONA(IEGLND(1),INVCN(1,IGLND)) .EQ.
     &                        ICONA(JJ,INVCN(IL,IGLND)))THEN
                              IEDGCT = IEDGCT + 1
                            END IF
                          ELSE
                            GO TO 380
                          END IF
 350                    CONTINUE
                        IF (IEDGCT .LT. 2)GO TO 380
 340                 CONTINUE
                     IF (NOWLND .EQ. IEGLND(1))THEN
                        IF ( KC .LT. 2) GO TO 380
                        NXTLND = IEGLND(2)
                        NXGLND = ICONA(NXTLND,INVCN(1,IGLND))
                     ELSE
                        NXTLND = IEGLND(1)
                        NXGLND = ICONA(NXTLND,INVCN(1,IGLND))
                     END IF

                     IF (INVLN(NXGLND) .LT. INVLN(IGLND))GO TO 380
                     IF (INVLN(NXGLND) .LT. 5)THEN
                        CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                       SOLEA,SOLENA,ITT,iblk)
                     ELSE
                        CALL EXTH(IGLND,INVCN,MAXLN,NXGLND,
     $                       INVLN(NXGLND),XA,YA,ZA,CNTRA,
     $                       SOLEA,SOLENA,ITT,iblk)
                     END IF
                     GO TO 200
                  END IF
 380              CONTINUE
                  IF (INVLN(IGLND) .LT. 5)THEN
                     CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                    SOLEA,SOLENA,ITT,iblk)
                  ELSE
                     CALL EXTH(IGLND,INVCN,MAXLN,IGLND,INVLN(IGLND),
     &                    XA,YA,ZA,CNTRA,SOLEA,SOLENA,ITT,iblk)
                  END IF
 200           CONTINUE

c*****

            ELSE IF (ITYPE .EQ. 6)THEN

c     Do for tet what you do for hex just not the same way

               DO 500 INOD = 1, NUMNDA
                  IGLND = NDLSTA(INOD)

c     First find elements connected to node - inverse connectivity
C     [NOTE: THIS DOES NOT SEEM TO BE USED...NOWLND?]
                  DO 510 I = 1, NNODES
                     IF (IGLND .EQ. ICONA(I,INVCN(1,IGLND)))THEN
                        NOWLND = I
                        GO TO 520
                     END IF
 510              CONTINUE
 520              CONTINUE

C     Less than 12 elements sharing IGLND, find the node of the
C     12 elements that connects with the maximum number of elements

C...Still not sure if this is correct for tets...
                  NDMAX = IGLND
                  IF (INVLN(IGLND) .LE. 12) THEN
                     DO 530 I = 1, INVLN(IGLND)
                        IEL = INVCN(I,IGLND)
                        DO 540 J = 1, NNODES
                           NOWNOD = ICONA(J,IEL)
                           IF (INVLN(NOWNOD) .GT. INVLN(NDMAX)) THEN
                              NDMAX = NOWNOD
                           END IF
 540                    CONTINUE
 530                 CONTINUE
                  ELSE
                     CALL EXTH(IGLND,INVCN,MAXLN,NDMAX,INVLN(NDMAX),
     &                    XA,YA,ZA,CNTRA,SOLEA,SOLENA,ITT,iblk)
                  END IF
                  IF (INVLN(NDMAX) .GT. 12 .and.
     $                 INVLN(IGLND) .GT.  4) THEN
                     CALL EXTH(IGLND,INVCN,MAXLN,IGLND,INVLN(IGLND),
     &                    XA,YA,ZA,CNTRA,SOLEA,SOLENA,ITT,iblk)
                  ELSE
                     CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),
     $                    SOLEA,SOLENA, ITT,iblk)
                  END IF
 500           CONTINUE
            END IF
            RETURN
            END
