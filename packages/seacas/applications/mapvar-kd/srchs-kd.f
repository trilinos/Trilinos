C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE SRCHS (
     *  NPSRF,   NFSRF,   LINKSRF, XYZSRF,
     *  NPTS,    XYZPTS,  TOLSRCH,
     *  NISR,    NRSR,    NISS,    NRSS,    ISRCHR,  RSRCHR,
     *  LIST,    IERR )

C-----------------------------------------------------------------------

C DESCRIPTION:

C THIS SUBROUTINE CALCULATES THE CLOSEST POINT PROBLEM
C BETWEEN NPTS POINTS AND NFSRF SURFACES AND RETURNS RESULTS OF
C SEARCH IN ISRCHR,RSRCHR

C USED HERE FOR FINDING LOCATION OF EITHER NODE OR ELEMENT CENTROID
C FROM MESH-B IN SHELL ELEMENT OF MESH-A

C-----------------------------------------------------------------------

C FORMAL PARAMETERS

C MEMORY      : P=PERMANENT, S=SCRATCH
C NAME        : IMPLICIT A-H,O-Z REAL, I-N INTEGER
C TYPE        : INPUT_STATUS/OUTPUT_STATUS (I=INPUT,O=OUTPUT,P=PASSED,
C               U=UNMODIFIED,-=UNDEFINED)
C DESCRIPTION : DESCRIPTION OF VARIABLE

C-----------------------------------------------------------------------

C CALLING ARGUMENTS:

C MEMORY NAME     TYPE   DESCRIPTION
C ---    ----     ---    -----------
C  P     NPSRF    I/U    NUMBER OF POINTS THAT DEFINE THE SURFACE
C  P     NFSRF    I/U    NUMBER OF SURFACES
C  P     LINKSRF  I/U    CONNECTIVITY OF SURFACES OF SIZE (4*NFSRF),
C                        NUMBERS REFER TO LOCATIONS IN XYZSRF ARRAY
C  P     XYZSRF   I/U    XYZ COORDS OF POINTS DEFINING SURFACE
C  P     NPTS     I/U    NUMBER OF POINTS TO BE SEARCHED
C  P     XYZPTS   I/U    XYZ COORDS OF POINTS TO BE SEARCHED
C  P     TOLSRCH  I/U    PROXIMITY TOLERANCE FOR POINT-TO-SURFACE SEARCH
C  P     NISR     I/U    NUMBER OF INTEGER SEARCH RESULTS (>=1)
C  P     NRSR     I/U    NUMBER OF REAL SEARCH RESULTS (>=4)
C  P     NISS     I/U    NUMBER OF INTEGER SEARCH SCRATCH (=5)
C  P     NRSS     I/U    NUMBER OF REAL SEARCH SCRATCH (=10)
C  P     ISRCHR   -/O    INTEGER SEARCH RESULTS
C  P     RSRCHR   -/O    REAL SEARCH RESULTS
C  S     LIST     -/-    LIST OF POTENTIAL CONTACTS FOR A SURFACE
C  S     IND      -/-    INDEX ARRAY INDICATING COORD. ORDER OF EACH POINT
C  S     IRNK     -/-    RANK ARRAY
C  S     IRNK2    -/-    RANK ARRAY (INDIRECT)
C  S     INDX     -/P    INTERMEDIATE LIST OF POTENTIAL PAIRS (CONSTRUCTED
C                          AFTER THE INTERSECTION OF TWO OF THREE LISTS)
C  S     ILO      -/-    SEARCH BOX MINIMUM INDEX
C  S     IUP      -/-    SEARCH BOX MAXIMUM INDEX
C  P     IERR     -/O    ERROR FLAG

C-----------------------------------------------------------------------

      include 'tapes.blk'
      include 'debg.blk'

C INPUT/OUTPUT ARRAYS

C ... Needed to interact with C routines on 64-bit systems which have
C     8-byte integers in Fortran and 4-byte integers in C.
      INTEGER*4 NPTS4, NDIM4, NLIST, LIST(*)

      DIMENSION
     *  LINKSRF(4,NFSRF),   XYZSRF(NPSRF,3),
     *  XYZPTS(NPTS,3),
     *  ISRCHR(NISR,NPTS),  RSRCHR(NRSR,NPTS)
      DIMENSION
     *  XMIN(3), XMAX(3), GXMIN(3), GXMAX(3)

C ISRCHR and RSRCHR must be initialized to zero

      DO 1 I = 1, NPTS
        DO 2 J = 1, NISR
          ISRCHR(J,I) = 0
 2      CONTINUE
        DO 3 K = 1, NRSR
          RSRCHR(K,I) = 0.
 3      CONTINUE
 1    CONTINUE

      IF( NISR .LT. 1 .OR. NRSR .LT. 4 .OR. NISS .LT. 5 .OR.
     *  NRSS .LT. 10 )THEN
        IERR = 1
        RETURN
      ENDIF
C DIMENSION OF COORDINATES
      NDIM = 3
      NDIM4 = 3
      NPTS4 = NPTS

c ... Calculate min/max extents of all points...
      GXMIN(1) = XYZPTS(1,1)
      GXMAX(1) = XYZPTS(1,1)
      GXMIN(2) = XYZPTS(1,2)
      GXMAX(2) = XYZPTS(1,2)
      GXMIN(3) = XYZPTS(1,3)
      GXMAX(3) = XYZPTS(1,3)
      DO 10 I=2,NPTS
        GXMIN(1) = MIN(GXMIN(1), XYZPTS(I,1))
        GXMIN(2) = MIN(GXMIN(2), XYZPTS(I,2))
        GXMIN(3) = MIN(GXMIN(3), XYZPTS(I,3))

        GXMAX(1) = MAX(GXMAX(1), XYZPTS(I,1))
        GXMAX(2) = MAX(GXMAX(2), XYZPTS(I,2))
        GXMAX(3) = MAX(GXMAX(3), XYZPTS(I,3))
 10   CONTINUE

C Build KD Tree
      if (idebug .ge. 2) then
        call excpus(t1)
        write(nout, *) '        In kdBuildTree ', npts4
        write(ntpout, *) '      In kdBuildTree', npts4
      end if

      call kdbuildtree(xyzpts, npts4, ndim4)

      if (idebug .ge. 2) then
        call excpus(t2)
        write(nout, *) '        Out of kdBuildTree', t2-t1
        write(ntpout, *) '      Out of kdBuildTree', t2-t1
      end if

C LOOP OVER SURFACES AND SEARCH FOR POINTS WITHIN CAPTURE BOX
      qt = 0.0
      j = 0
      nlistt = 0
      nlistp = 0
      nskip  = 0
      call excpus(t3p)

C LOOP OVER SURFACES AND SEARCH FOR POINTS WITHIN CAPTURE BOX
      DO 100 IFSRF = 1, NFSRF
        NI = LINKSRF(1,IFSRF)
        NJ = LINKSRF(2,IFSRF)
        NK = LINKSRF(3,IFSRF)
        NL = LINKSRF(4,IFSRF)

        XMINMS = MIN(XYZSRF(NI,1),XYZSRF(NJ,1),
     *               XYZSRF(NK,1),XYZSRF(NL,1))
        XMAXMS = MAX(XYZSRF(NI,1),XYZSRF(NJ,1),
     *               XYZSRF(NK,1),XYZSRF(NL,1))
        YMINMS = MIN(XYZSRF(NI,2),XYZSRF(NJ,2),
     *               XYZSRF(NK,2),XYZSRF(NL,2))
        YMAXMS = MAX(XYZSRF(NI,2),XYZSRF(NJ,2),
     *               XYZSRF(NK,2),XYZSRF(NL,2))
        ZMINMS = MIN(XYZSRF(NI,3),XYZSRF(NJ,3),
     *               XYZSRF(NK,3),XYZSRF(NL,3))
        ZMAXMS = MAX(XYZSRF(NI,3),XYZSRF(NJ,3),
     *               XYZSRF(NK,3),XYZSRF(NL,3))
        XDEL = XMAXMS - XMINMS
        YDEL = YMAXMS - YMINMS
        ZDEL = ZMAXMS - ZMINMS
        TOLER = TOLSRCH * MAX(XDEL,YDEL,ZDEL)
        XMIN(1) = XMINMS - TOLER
        XMAX(1) = XMAXMS + TOLER
        XMIN(2) = YMINMS - TOLER
        XMAX(2) = YMAXMS + TOLER
        XMIN(3) = ZMINMS - TOLER
        XMAX(3) = ZMAXMS + TOLER

C ... Build a list of points in the query region...
        j = j + 1
C ... Skip past points that are outside search domain.
C     This is much faster than the failed search in kdrectquery, but
C     does add overhead if all points are in the domain....

        if (xmax(1) .lt. gxmin(1) .or. xmin(1) .gt. gxmax(1) .or.
     *      xmax(2) .lt. gxmin(2) .or. xmin(2) .gt. gxmax(2) .or.
     *    xmax(3) .lt. gxmin(3) .or. xmin(3) .gt. gxmax(3)) then
          nskip = nskip + 1
          nlist = 0
        else
          call excpus(t4)
          call kdrectquery(xyzpts,npts4,ndim4, xmin,xmax, list,nlist)
          call excpus(t5)
          qt = qt + (t5-t4)
        end if

C ... Debugging statistics.  Also give user an idea that calculation is progressing...
        if (idebug .ge. 2) then
          nlistp = nlistp + nlist
          nlistt = nlistt + nlist
          if (j .eq. 100000 .or. ifsrf .eq. nfsrf) then
            call excpus(t3)
            write (nout, 190) ifsrf, nfsrf, t3-t3p, t3-t2,
     *        j/(t3-t3p), ifsrf/(t3-t2), nlistp, nlistt
            write (ntpout, 190) ifsrf, nfsrf, t3-t3p, t3-t2,
     *        j/(t3-t3p), ifsrf/(t3-t2), nlistp, nlistt
            j = 0
            t3p = t3
            nlistp = 0
          end if
        end if

        DO 140 K = 1, NLIST
          lval = list(k)
          CALL SHLSRC(
     *      NDIM,     NPTS,     NPSRF,    NFSRF,    NISR,
     *      NRSR,     NRSS,     XYZSRF,   XYZPTS,   LINKSRF,
     *      ISRCHR,   RSRCHR,   LVAL,     IFSRF,    TOLSRCH,  IERR   )

 140    CONTINUE
 100  CONTINUE

C ... More debugging stats
      if (idebug .ge. 2) then
        call excpus(t3)
        write(nout, *) '  Finish Queries', t3-t2, NFSRF/(t3-t2), qt
        write(ntpout, *) '  Finish Queries', t3-t2, NFSRF/(t3-t2), qt
        write(nout, *) '  Matches = ',nlistt,' Rate = ',nlistt/(t3-t2)
        write(ntpout, *) '  Matches = ',nlistt,' Rate = ',nlistt/(t3-t2)
      end if

      call kdkilltree()

 190  format(1x,i9,'/',i9,' T= ',1pe10.3,'/',1pe10.3,
     *  ', R= ',1pe10.3,'/',1pe10.3,' M = ',i10,'/',i10)

      RETURN
      END

