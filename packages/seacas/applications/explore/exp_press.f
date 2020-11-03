C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRESS (OPTION, NOUT, NUMESS, LISESS, LESSEL, LESSNL,
     &  IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTSESS, FACESS,
     *  NAME,  nvar, namvar, isvok, lisvar, NDFSID, NODSID,
     *  MAPEL, MAPND, DOMAPE, DOMAPN)
C=======================================================================

C     --*** PRESS *** (EXPLORE) Display database element side set
C     --
C     --PRESS displays the element side sets.
C     --
C     --Parameters:
C     --   OPTION - IN - '*' to print all, else print options:
C     --      ' ' - set summary (number of nodes, etc)
C     --      'E' - elements in set (may not be combined with 'N' or 'F')
C     --      'N' - nodes in set (may not be combined with 'E')
C     --      'F' - distribution factors for set (may not be combined with 'E')
C     --   NOUT - IN - the output file, <=0 for standard
C     --   NUMESS - IN - the number of element side sets
C     --   LISESS - IN - the indices of the selected element side sets
C     --   LESSEL - IN - the number of elements for all sets
C     --   LESSNL - IN - the number of nodes for all sets
C     --   IDESS - IN - the element side set ID for each set
C     --   NEESS - IN - the number of elements for each set
C     --   NNESS - IN - the number of nodes for each set
C     --   IXEESS - IN - the index of the first element for each set
C     --   IXNESS - IN - the index of the first node for each set
C     --   LTEESS - IN - the elements for all sets
C     --   LTSESS - IN - the element sides for all sets
C     --   FACESS - IN - the distribution factors for all sets
C     --   NVAR - IN - the number of variables
C     --   NAMVAR - IN - the names of the variables
C     --   NDFSID - IN - the number of df per face
C     --   ISVOK  - IN - the variable truth table;
C     --      variable i of set j exists iff ISVOK(i,j) is NOT 0
C     --   LISVAR  - SCRATCH - size = NVAR (if 'V' in OPTION)

      include 'exodusII.inc'
      INCLUDE 'exp_dbase.blk'

      CHARACTER*(*) OPTION
      INTEGER LISESS(0:*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      INTEGER NDFSID(*)
      INTEGER NODSID(*)
      REAL FACESS(*)
      CHARACTER*(*) NAME(*)
      CHARACTER*(*) NAMVAR(*)
      INTEGER ISVOK(NVAR,*)
      INTEGER LISVAR(*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)
      LOGICAL DOMAPE, DOMAPN

      LOGICAL ALLSAM
      LOGICAL DOELE, DONOD, DOFAC, DOVTBL
      CHARACTER*20 STRA, STRB, STRC

      INTEGER GETPRC, PRTLEN
      CHARACTER*128 FMT1, FMTE, FMTM

      INTEGER NODES(27)
      REAL    FACTORS(27)

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMTE, 10055) FMT1(:LFMT)
      WRITE(FMTM, 10100) FMT1(:LFMT)

      DOELE  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0))
      DONOD  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOFAC  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))

      IF (NOUT .GT. 0) THEN
        IF (DOELE) THEN
          WRITE (NOUT, 10020) 'ELEMENT LIST'
        ELSE IF (DONOD .AND. DOFAC) THEN
          WRITE (NOUT, 10020) 'ELEM.FACE LIST AND DISTRIBUTION FACTORS'
        ELSE IF (DONOD) THEN
          WRITE (NOUT, 10020) 'ELEM.FACE LIST'
        ELSE IF (DOFAC) THEN
          WRITE (NOUT, 10020) 'DISTRIBUTION FACTORS'
        ELSE
          WRITE (NOUT, 10020)
        END IF
        if (domape) then
          write (nout, 10025)
        end if
      ELSE
        if (domape) then
          write (*, 10025)
        end if
      END IF

      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, *)
      ELSE
        WRITE (*, *)
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NUMESS
10000 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LESSEL
10010 FORMAT ('(index=', I12, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)
      WRITE (STRC, 10010, IOSTAT=IDUM) LESSNL
      CALL PCKSTR (1, STRC)
      LSTRC = LENSTR (STRC)

      DO 100 IX = 1, LISESS(0)
        IESS = LISESS(IX)
        WRITE (STRA, 10000, IOSTAT=IDUM) IESS
        CALL PCKSTR (1, STRA)
        WRITE (STRB, 10010, IOSTAT=IDUM) IXEESS(IESS)
        CALL PCKSTR (1, STRB)
        WRITE (STRC, 10010, IOSTAT=IDUM) IXNESS(IESS)
        CALL PCKSTR (1, STRC)
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10030, IOSTAT=IDUM)
     &      IDESS(IESS), STRA(:LSTRA),
     &      NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $      NAME(IESS)(:LENSTR(NAME(IESS)))
        ELSE
          WRITE (*, 10030, IOSTAT=IDUM)
     &      IDESS(IESS), STRA(:LSTRA),
     &      NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $      NAME(IESS)(:LENSTR(NAME(IESS)))
        END IF

        IF (DOELE .AND. (NEESS(IESS) .GT. 0)) THEN
          IS = IXEESS(IESS)
          IE = IS + NEESS(IESS) - 1
          IF (NOUT .GT. 0) THEN
            if (domape) then
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            else
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (LTEESS(I), I=IS,IE)
            end if
          ELSE
            if (domape) then
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            else
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (LTEESS(I), I=IS,IE)
            end if
          END IF
        END IF

C     ... This used to print the nodes of the sideset faces, it now prints
C     the local element faces of the sideset
        IF (DONOD .AND. (NEESS(IESS) .GT. 0)) THEN
          IS = IXEESS(IESS)
          IE = IS + NEESS(IESS) - 1
          IF (NOUT .GT. 0) THEN
            if (domape) then
              WRITE (NOUT, 10045, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)),LTSESS(I), I=IS,IE)
            else
              WRITE (NOUT, 10045, IOSTAT=IDUM)
     &          (LTEESS(I),LTSESS(I), I=IS,IE)
            endif
          ELSE
            if (domape) then
              WRITE (*, 10045, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)),LTSESS(I), I=IS,IE)
            else
              WRITE (*, 10045, IOSTAT=IDUM)
     &          (LTEESS(I),LTSESS(I), I=IS,IE)
            end if
          END IF
        END IF

        IF (DOVTBL) THEN
          NSEL = 0
          mxnam = 0
          DO 30 I = 1, NVAR
            IF (ISVOK(I,IESS) .NE. 0) THEN
              NSEL = NSEL + 1
              LISVAR(NSEL) = I
              LNAM = lenstr(namvar(i))
              mxnam = max(lnam, mxnam)
            END IF
 30       CONTINUE

          if (nsel .gt. 0) then
            IF (NOUT .GT. 0) THEN
              if (mxnam .gt. 24) then
                WRITE (NOUT, 10090, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              else
                WRITE (NOUT, 10070, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              endif
              write (nout, 10080)
            ELSE
              if (mxnam .gt. 24) then
                WRITE (*, 10090, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              else
                WRITE (*, 10070, IOSTAT=IDUM)
     &            (NAMVAR(LISVAR(I))(:mxnam), I=1,NSEL)
              endif
              write (*, 10080)
            END IF
          end if
        END IF

        IF (DOFAC) THEN
          IF (NNESS(IESS) .GT. 0) THEN
            IS = IXNESS(IESS)
            IE = IS + NNESS(IESS) - 1
C     ... See if all values are the same
            val = facess(is)
            allsam = .TRUE.
            do 50 i=is+1, ie
              if (facess(i) .ne. val) then
                allsam = .FALSE.
                go to 90
              end if
 50         continue
 90         continue
            if (allsam) then
              IF (NOUT .GT. 0) THEN
                WRITE (NOUT, FMTE, IOSTAT=IDUM) VAL
              ELSE
                WRITE (*, FMTE, IOSTAT=IDUM) VAL
              END IF
            else
C ... Get the number of df/nodes per face...
C ... NOTE: facess is contiguous over all sidesets,
C           nodsid is only for the current sideset.
              call exgssn(ndb, idess(iess), ndfsid, nodsid, ierr)
              ISE = IXEESS(IESS)
              IEE = ISE + NEESS(IESS) - 1
              IDS = IS
              idf = 1
              idn = 1
              do i=ise, iee
                NDFPE = ndfsid(idf)
                idf=idf+1
                if (domape) then
                  iel = mapel(lteess(i))
                else
                  iel = lteess(i)
                end if
                do j=1,ndfpe
                  factors(j) = facess(j+ids-1)
                  nodes(j) = nodsid(j+idn-1)
                end do

                if (domapn) then
                  do j=1,ndfpe
                    nodes(j) = mapnd(nodes(j))
                  end do
                end if

                if (nout .gt. 0) then
                  write (nout, FMTM) iel, ltsess(i),
     *              (nodes(j), factors(j),j=1,ndfpe)
                else
                  write (*, FMTM) iel, ltsess(i),
     *              (nodes(j), factors(j),j=1,ndfpe)
                end if
                ids = ids + ndfpe
                idn = idn + ndfpe
              end do
            end if
          else
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10060, IOSTAT=IDUM)
            ELSE
              WRITE (*, 10060, IOSTAT=IDUM)
            END IF
          end if
        END IF
 100  CONTINUE

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)
10020 FORMAT (/, 1X, 'ELEMENT SIDE SETS', :, ' - ', A)
10025 FORMAT (1x, 'Element Ids are Global')
10030 FORMAT (1X, 'Set', I12, 1X, A, ':',
     &  I12, ' elements', 1X, A,
     &  I12, ' nodes/df', 1X, A, ' name = "',A,'"')
10040 FORMAT ((1X, 8I12))
10045 FORMAT ((1X, 8(I12,'.',I1)))
10055 FORMAT ('(10x, ''All distribution factors are equal to ''', A,')')
10060 FORMAT (10x, 'Distribution factors not stored in file.')
10070 FORMAT ((2x,4(2X, A)))
10090 FORMAT ((2x,3(2X, A)))
10080 FORMAT (1X)
10100 FORMAT ('((1X, I12,''.'',I1,8x,(8(1x,I12,1x,',A,')))))')
      END
