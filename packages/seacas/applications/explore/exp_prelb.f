C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRELB (OPTION, NOUT, NELBLK, NLISEL, LISEL,
     &  IDELB, LENE, NUMLNK, NUMATR, LINK, ATRIB,
     &  EBTYPE, EBNAME, NVAREL, NAMEEV, ISEVOK, LISEV,
     *  ATNAME, MAPN, DOMAPN, MAPE, DOMAPE)
C=======================================================================

C   --*** PRELB *** (BLOT) Display database element blocks
C   --
C   --PRELB displays the element blocks.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' to print block summary
C   --      'N' to print element block name
C   --      'V' to print element variable truth table
C   --      'C' to print connectivity
C   --      'A' to print attributes
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NELBLK - IN - the number of element blocks
C   --   NLISEL - IN - the number of selected elements by block
C   --   LISEL - IN - the indices of the selected elements by block
C   --   IDELB - IN - the element block ID for each block
C   --   LENE - IN - the cumulative element counts by element block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the connectivity array for all blocks
C   --   ATRIB - IN - the attribute array for all blocks
C   --   EBTYPE - IN - the names of the element block types
C   --   NVAREL - IN - the number of element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   LISEV - SCRATCH - size = NVAREL (if 'V' in OPTION)

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      INTEGER NLISEL(0:NELBLK)
      INTEGER LISEL(0:*)
      INTEGER IDELB(*)
      INTEGER LENE(0:NELBLK)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(*) EBNAME(*)
      CHARACTER*(*) NAMEEV(*)
      CHARACTER*(*) ATNAME(*)
      INTEGER ISEVOK(NVAREL,*)
      INTEGER LISEV(*)
      INTEGER MAPN(*), MAPE(*)
      LOGICAL DOMAPN, DOMAPE
      LOGICAL ISABRT
      LOGICAL DONAM, DOVTBL, DOCONN, DOATR
      LOGICAL BLK1
      CHARACTER*40 STRA, STRB

      DONAM  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))
      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))

      DONAM = .TRUE.
      BLK1 = .TRUE.
      IF (NOUT .GT. 0) THEN
        IF (DOCONN .AND. DOATR) THEN
          WRITE (NOUT, 10020) 'CONNECTIVITY and ATTRIBUTES'
        ELSE IF (DOCONN) THEN
          WRITE (NOUT, 10020) 'CONNECTIVITY'
        ELSE IF (DOATR) THEN
          WRITE (NOUT, 10020) 'ATTRIBUTES'
        ELSE
          WRITE (NOUT, 10020)
        END IF
      END IF

      if (nout .gt. 0) then
        if (domape) write (nout, 10060) 'Element'
        if (domapn) write (nout, 10060) 'Nodal'
      else
        if (domape) write (*, 10060) 'Element'
        if (domapn) write (*, 10060) 'Nodal'
      end if

      WRITE (STRA, 10000, IOSTAT=IDUM) NELBLK
10000 FORMAT ('(#', I12, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LENE(NELBLK), LENE(NELBLK)
10010 FORMAT ('(', I12, '..', I12, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      ISLNK = 1
      ISATR = 1
      ISATN = 1

      DO 110 IELB = 1, NELBLK
        NUME = LENE(IELB) - LENE(IELB-1)
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
          numnod = numlnk(ielb)
        else
          numnod = numlnk(ielb)*nume
        end if
        IF (ISABRT ()) RETURN
        IF (NLISEL(IELB) .GT. 0) THEN
          IEL = LENE(IELB-1)+1
          LEL = LENE(IELB)

          IF (BLK1 .OR. DOCONN .OR. DOATR) THEN
            BLK1 = .FALSE.
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, *)
            ELSE
              WRITE (*, *)
            END IF
          END IF

          WRITE (STRA, 10000, IOSTAT=IDUM) IELB
          CALL PCKSTR (1, STRA)
          WRITE (STRB, 10010, IOSTAT=IDUM) IEL, LEL
          CALL PCKSTR (1, STRB)
          IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &        IDELB(IELB), STRA(:LSTRA),
     &        NUME, STRB(:LSTRB), NUMLNK(IELB), NUMATR(IELB)
          ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &        IDELB(IELB), STRA(:LSTRA),
     &        NUME, STRB(:LSTRB), NUMLNK(IELB), NUMATR(IELB)
          END IF

          IF (DONAM) THEN
            LNAM = LENSTR(EBTYPE(IELB))
            LNM =  LENSTR(EBNAME(IELB))
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10040) EBNAME(IELB)(:LNM),
     $          EBTYPE(IELB)(:LNAM)
            ELSE
              WRITE (*, 10040) EBNAME(IELB)(:LNM), EBTYPE(IELB)(:LNAM)
            END IF
            if (numatr(ielb) .gt. 0) then
              mxnam = 0
              DO I = 1, numatr(ielb)
                LNAM = lenstr(atname(isatn+i-1))
                mxnam = max(lnam, mxnam)
              END DO
              if (mxnam .gt. 1) then
                IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10090, IOSTAT=IDUM)
     &              (ATNAME(I)(:mxnam),
     *              I=ISATN,ISATN+NUMATR(IELB)-1)
                ELSE
                  WRITE (*, 10090, IOSTAT=IDUM)
     &              (ATNAME(I)(:mxnam),
     *              I=ISATN,ISATN+NUMATR(IELB)-1)
                END IF
              end if
            end if
          END IF

          IF (DOVTBL) THEN
            NSEL = 0
            mxnam = 0
            DO 100 I = 1, NVAREL
              IF (ISEVOK(I,IELB) .NE. 0) THEN
                NSEL = NSEL + 1
                LISEV(NSEL) = I
                LNAM = lenstr(nameev(i))
                mxnam = max(lnam, mxnam)
              END IF
 100        CONTINUE

            IF (NOUT .GT. 0) THEN
              if (mxnam .gt. 24) then
                WRITE (NOUT, 10070, IOSTAT=IDUM)
     &            (NAMEEV(LISEV(I))(:mxnam), I=1,NSEL)
              else
                WRITE (NOUT, 10050, IOSTAT=IDUM)
     &            (NAMEEV(LISEV(I))(:mxnam), I=1,NSEL)
              endif
              write (nout, 10080)
            ELSE
              if (mxnam .gt. 24) then
                WRITE (*, 10070, IOSTAT=IDUM)
     &            (NAMEEV(LISEV(I))(:mxnam), I=1,NSEL)
              else
                WRITE (*, 10050, IOSTAT=IDUM)
     &            (NAMEEV(LISEV(I))(:mxnam), I=1,NSEL)
              endif
              write (*, 10080)
            END IF
          END IF

          IF (DOCONN .OR. DOATR) THEN
            if (ebtype(ielb) .eq. 'nsided' .or.
     *        ebtype(ielb) .eq. 'NSIDED') THEN
              CALL PREBN (OPTION, NOUT, IELB,
     *          nume, IDELB(IELB),
     *          IEL-1, NLISEL(IELB), LISEL(IEL),
     *          NUMATR(IELB),
     &          LINK(ISLNK), ATRIB(ISATR), max(1,numatr(ielb)),
     *          MAPN, DOMAPN, MAPE, DOMAPE)
            else
              CALL PREB1 (OPTION, NOUT, IEL-1,
     &          NLISEL(IELB), LISEL(IEL), NUMLNK(IELB), NUMATR(IELB),
     &          LINK(ISLNK), ATRIB(ISATR), max(1,numatr(ielb)),
     *          MAPN, DOMAPN, MAPE, DOMAPE)
            end if
            IF (ISABRT ()) RETURN
          END IF
        END IF

        IF (DOCONN) ISLNK = ISLNK + numnod
        IF (DOATR)  ISATR = ISATR + NUME * NUMATR(IELB)
        ISATN = ISATN + NUMATR(IELB)
 110  CONTINUE

      RETURN

10020 FORMAT (/, 1X, 'ELEMENT BLOCKS', :, ' - ', A)
C ... NOTE: Normal elements have a a mximum of <100 nodes, but superelements
C           Or other "strange" elements may have lots of nodes; keep format high
10030 FORMAT (1X, 'Block', I12, 1X, A, ':',
     &  I12, ' elements', 1X, A,
     &  I6, '-node', I4, ' attributes')
10040 FORMAT (2X, 'Element block name = "',A
     $  ,'", Element type = "', A, '"')
10050 FORMAT ((2x,4(2X, A)))
10060 FORMAT (1x, A, ' IDs are global')
10070 FORMAT ((2x,3(2X, A)))
10080 FORMAT (1X)
10090 FORMAT (2X, 'Attributes: ', 10(2X, A))
      END
