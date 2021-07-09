C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C   met driver is the file device
C   x11 driver is the terminal device
      SUBROUTINE VDESCP(ESCPCD,N,ARGS)
      INTEGER ESCPCD,SUPPRT,TRMSUP,FILSUP
      REAL ARGS(*)
      LOGICAL TERMON,FILEON,TERMUP,FILEUP

      real locx, locy
      real clrary(*)
      integer clrmod
      integer index, indexa(*)
      real attarr(*)
      real aspect
      integer justif
      real linwth
      integer chars(*), color
      real xarray(*), yarray(*)
      real inten

#ifdef Build64
      integer*4 escpcd4, n4
      real*4 args4(100)
      real*4 aspect4
      integer*4 justif4
      integer*4 trmsup4, filsup4, supprt4, esc4
      integer*4 istat4
      integer*4 char4
      real*4 x4, y4
      integer*4 num4
      integer*4 index4, index4a(256), clrmod4
      real*4 clrary4(768)
      real*4 value4
      real*4 xndc4, yndc4
      real*4 attarr4(7)
      integer*4 color4
      real*4 inten4
      real*4 locx4, locy4
      integer*4 linsty4
      real*4 linwth4
      real*4 xa4(1024), ya4(1024)
      integer*4 npts4
      integer*4 length4
      integer*4 chars4(136)
#endif

      DATA TERMON,FILEON,TERMUP,FILEUP /2*.TRUE.,2*.TRUE./

      IF( ESCPCD.GE.10000 .AND. ESCPCD.LE.10003 ) THEN

C Special Escape Codes for Dual Device Control:
        KSTAT = ESCPCD - 10000
        TERMON = (KSTAT.EQ.1 .OR. KSTAT.EQ.3) .AND. TERMUP
        FILEON = (KSTAT.EQ.2 .OR. KSTAT.EQ.3) .AND. FILEUP
      ELSE
        IF(TERMON) THEN
#ifdef Build64
C ... For X11 driver, there is only a single arg ever used and it
C     is write only...
          escpcd4 = escpcd
          n4 = n
          CALL wx11es(ESCPCD4,N4,ARGS4)
          if (escpcd .eq. 3501) args(1) = args4(1)
#else
          CALL wx11es(ESCPCD,N,ARGS)
#endif
        END IF
        IF(FILEON) THEN
C ... This probably needs to be fixed....
          CALL wmetes(ESCPCD,N,ARGS)
        END IF
      END IF
      RETURN

      ENTRY VDIQES(ESCPCD,SUPPRT)

C Special Escape Codes for Dual Device Control:
      IF( ESCPCD.EQ.10000 ) THEN
        SUPPRT = 1
      ELSE IF ( ESCPCD.EQ.10001 ) THEN
#ifdef Build64
        esc4 = 10000
        CALL wx11ie(esc4,TRMSUP4)
        TRMSUP = TRMSUP4
#else
        CALL wx11ie(10000,TRMSUP)
#endif
        SUPPRT = 1 - TRMSUP
      ELSE IF ( ESCPCD.EQ.10002 ) THEN
#ifdef Build64
        esc4 = 10000
        CALL wmetie(esc4,FILSUP4)
        FILSUP = FILSUP4
#else
        CALL wmetie(10000,FILSUP)
#endif
        SUPPRT = 1 - FILSUP
      ELSE IF ( ESCPCD.EQ.10003 ) THEN
#ifdef Build64
        esc4 = 10000
        CALL wx11ie(esc4,TRMSUP4)
        TRMSUP = TRMSUP4
#else
        CALL wx11ie(10000,TRMSUP)
#endif
#ifdef Build64
        esc4 = 10000
        CALL wmetie(esc4,FILSUP4)
        FILSUP = FILSUP4
#else
        CALL wmetie(10000,FILSUP)
#endif
        SUPPRT = (1 - TRMSUP) * (1 - FILSUP)
      ELSE
        IF(FILEON) THEN
#ifdef Build64
          escpcd4 = escpcd
          CALL wmetie(ESCPCD4,SUPPRT4)
          SUPPRT = SUPPRT4
#else
          CALL wmetie(ESCPCD,SUPPRT)
#endif
        END IF
        IF(TERMON) THEN
#ifdef Build64
          escpcd4 = escpcd
          CALL wx11ie(ESCPCD4,SUPPRT4)
          SUPPRT = SUPPRT4
#else
          CALL wx11ie(ESCPCD,SUPPRT)
#endif
        END IF
      END IF
      RETURN

      ENTRY VDINIT(ASPECT,JUSTIF)

C Check if terminal and hardcopy devices are null:
#ifdef Build64
      ESC4 = 10000
      CALL wx11ie( ESC4,ISTAT4 )
      ISTAT = ISTAT4
#else
      CALL wx11ie( 10000,ISTAT )
#endif
      TERMUP = ISTAT .EQ. 0
      TERMON = TERMON .AND. TERMUP
#ifdef Build64
      ESC4 = 10000
      CALL wmetie( ESC4,ISTAT4 )
      ISTAT = ISTAT4
#else
      CALL wmetie( 10000,ISTAT )
#endif
      FILEUP = ISTAT .EQ. 0
      FILEON = FILEON .AND. FILEUP

C Stop program if both devices are null.
      IF(.NOT.(TERMUP.OR.FILEUP)) STOP 'No graphics device assigned.'

C      print *, 'About to init x11'
      IF(TERMON) THEN
#ifdef Build64
        ASPECT4 = ASPECT
        JUSTIF4 = JUSTIF
        CALL wx11nt(ASPECT4,JUSTIF4)
#else
        CALL wx11nt(ASPECT,JUSTIF)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        ASPECT4 = ASPECT
        JUSTIF4 = JUSTIF
c        print *, 'About to init metafile 64 ', aspect4, justif4
        CALL wmetnt(ASPECT4,JUSTIF4)
#else
c        print *, 'About to init metafile', aspect, justif
        CALL wmetnt(ASPECT,JUSTIF)
#endif
c        print *, 'Done with init', FILEON
      END IF
      RETURN

      ENTRY VDAABU(BTNNUM)
      IF(FILEON) CALL wmetbu(BTNNUM)
      IF(TERMON) THEN
#ifdef Build64
C ... Not implemented in X11 driver
        CALL wx11bu(BTNNUM)
#else
        CALL wx11bu(BTNNUM)
#endif
      END IF
      RETURN

      ENTRY VDABGL(BTNNUM,X,Y)
      IF(FILEON) CALL wmetbl(BTNNUM,X,Y)
      IF(TERMON) THEN
#ifdef Build64
C ... Not implemented in X11 driver
        CALL wx11bl(BTNNUM,X,Y)
#else
        CALL wx11bl(BTNNUM,X,Y)
#endif
      END IF
      RETURN

      ENTRY VDAKGL(CHAR,X,Y)
      IF(FILEON) THEN
#ifdef Build64
        CALL wmetkl(CHAR4,X4,Y4)
        char = char4
        x = x4
        y = y4
#else
        CALL wmetkl(CHAR,X,Y)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        CALL wx11kl(CHAR4,X4,Y4)
        char = char4
        x = x4
        y = y4
#else
        CALL wx11kl(CHAR,X,Y)
#endif
      END IF
      RETURN

      ENTRY VDALOC(X,Y)
      IF(FILEON) THEN
#ifdef Build64
        CALL wmetlo(X4,Y4)
        x = x4
        y = y4
#else
        CALL wmetlo(X,Y)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        CALL wx11lo(X4,Y4)
        x = x4
        y = y4
#else
        CALL wx11lo(X,Y)
#endif
      END IF
      RETURN

      ENTRY VDBELL
      IF(TERMON) CALL wx11be
      IF(FILEON) CALL wmetbe
      RETURN

      ENTRY VDBUFL
      IF(TERMON) CALL wx11fl
      IF(FILEON) CALL wmetfl
      RETURN

      ENTRY VDFRAM(ITYPE)
C ... Doesn't do anything in X11 driver
      IF(TERMON) CALL wx11fr(ITYPE)
      IF(FILEON) CALL wmetfr(ITYPE)
      RETURN

      ENTRY VDIQCO(NUM,INDEXA,CLRARY,CLRMOD)
      IF(FILEON) THEN
#ifdef Build64
        num4 = num
        clrmod4 = clrmod
        do 20 i=1, num
          index4a(i) = indexa(i)
 20     continue
        CALL wmetic(NUM4,INDEX4a,CLRARY4,CLRMOD4)
        do 30 i=1, 3*num
          clrary(i) = clrary4(i)
 30     continue
#else
        CALL wmetic(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        num4 = num
        clrmod4 = clrmod
        do 21 i=1, num
          index4a(i) = indexa(i)
 21     continue
        CALL wx11ic(NUM4,INDEX4a,CLRARY4,CLRMOD4)
        do 31 i=1, 3*num
          clrary(i) = clrary4(i)
 31     continue
#else
        CALL wx11ic(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      RETURN

      ENTRY VDIQCP(X,Y)
      IF(FILEON) THEN
#ifdef Build64
        CALL wmetcp(X4,Y4)
        x = x4
        y = y4
#else
        CALL wmetcp(X,Y)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        CALL wx11cp(X4,Y4)
        x = x4
        y = y4
#else
        CALL wx11cp(X,Y)
#endif
      END IF
      RETURN

      ENTRY VDIQDC(INDEX,VALUE)
      IF(FILEON) THEN
#ifdef Build64
        index4 = index
        CALL wmetdc(INDEX4,VALUE4)
        value = value4
#else
        CALL wmetdc(INDEX,VALUE)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        index4 = index
        CALL wx11dc(INDEX4,VALUE4)
        value = value4
#else
        CALL wx11dc(INDEX,VALUE)
#endif
      END IF
      RETURN

      ENTRY VDIQND(XNDC,YNDC)
      IF(FILEON) THEN
#ifdef Build64
        CALL wmetid(XNDC4,YNDC4)
        xndc = xndc4
        yndc = yndc4
#else
        CALL wmetid(XNDC,YNDC)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        CALL wx11id(XNDC4,YNDC4)
        xndc = xndc4
        yndc = yndc4
#else
        CALL wx11id(XNDC,YNDC)
#endif
      END IF
      RETURN

      ENTRY VDIQOS(ATTARR)
      IF(FILEON) THEN
#ifdef Build64
        CALL wmetio(ATTARR4)
        do 41 i=1,7
          attarr(i) = attarr4(i)
 41     continue
#else
        CALL wmetio(ATTARR)
#endif
      END IF
      IF(TERMON) THEN
#ifdef Build64
        CALL wx11io(ATTARR4)
        do 40 i=1,7
          attarr(i) = attarr4(i)
 40     continue
#else
        CALL wx11io(ATTARR)
#endif
      END IF
      RETURN

      ENTRY VDLOGE(ERRNUM,ERRSEV)
C ... Doesn't do anything in X11 driver
      IF(TERMON) CALL wx11le(ERRNUM,ERRSEV)
      IF(FILEON) CALL wmetle(ERRNUM,ERRSEV)
      RETURN

      ENTRY VDSTBC(COLOR)
      IF(TERMON) THEN
#ifdef Build64
        color4 = color
        CALL wx11bc(COLOR4)
#else
        CALL wx11bc(COLOR)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        color4 = color
        CALL wmetbc(COLOR4)
#else
        CALL wmetbc(COLOR)
#endif
      END IF
      RETURN

      ENTRY VDSTCO(NUM,INDEXA,CLRARY,CLRMOD)
      IF(TERMON) THEN
#ifdef Build64
        num4 = num
        clrmod4 = clrmod
        j = 1
        do 50 i=1, num
          index4a(i) = indexa(i)
 50     continue
        do 55 i=1, num*3
          clrary4(i) = clrary(i)
 55     continue
        CALL wx11co(NUM4,INDEX4A,CLRARY4,CLRMOD4)
#else
        CALL wx11co(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        num4 = num
        clrmod4 = clrmod
        j = 1
        do 51 i=1, num
          index4a(i) = indexa(i)
 51     continue
        do 56 i=1, num*3
          clrary4(i) = clrary(i)
 56     continue
        CALL wmetco(NUM4,INDEX4A,CLRARY4,CLRMOD4)
#else
        CALL wmetco(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      RETURN

      ENTRY VDSTCS(YSIZE)
      IF (TERMON) THEN
#ifdef Build64
        Y4 = YSIZE
        CALL wx11cs(Y4)
#else
        CALL wx11cs(YSIZE)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        Y4 = YSIZE
        CALL wmetcs(Y4)
#else
        CALL wmetcs(YSIZE)
#endif
      END IF
      RETURN

      ENTRY VDSTFC(COLOR)
      IF (TERMON) THEN
#ifdef Build64
        color4 = color
        CALL wx11fc(COLOR4)
#else
        CALL wx11fc(COLOR)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        color4 = color
        CALL wmetfc(COLOR4)
#else
        CALL wmetfc(COLOR)
#endif
      END IF
      RETURN

      ENTRY VDSTIN(INTEN)
      IF (TERMON) THEN
#ifdef Build64
        inten4 = inten
        CALL wx11in(INTEN4)
#else
        CALL wx11in(INTEN)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        inten4 = inten
        CALL wmetin(INTEN4)
#else
        CALL wmetin(INTEN)
#endif
      END IF
      RETURN

      ENTRY VDSTLA(LOCX,LOCY)
      IF (TERMON) THEN
#ifdef Build64
        locx4 = locx
        locy4 = locy
        CALL wx11la(LOCX4,LOCY4)
#else
        CALL wx11la(LOCX,LOCY)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        locx4 = locx
        locy4 = locy
        CALL wmetla(LOCX4,LOCY4)
#else
        CALL wmetla(LOCX,LOCY)
#endif
      END IF
      RETURN

      ENTRY VDSTLS(LINSTY)
      IF (TERMON) THEN
#ifdef Build64
        linsty4 = linsty
        CALL wx11ls(LINSTY4)
#else
        CALL wx11ls(LINSTY)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        linsty4 = linsty
        CALL wmetls(LINSTY4)
#else
        CALL wmetls(LINSTY)
#endif
      END IF
      RETURN

      ENTRY VDSTLW(LINWTH)
      IF (TERMON) THEN
#ifdef Build64
        linwth4 = linwth
        CALL wx11lw(LINWTH4)
#else
        CALL wx11lw(LINWTH)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        linwth4 = linwth
        CALL wmetlw(LINWTH4)
#else
        CALL wmetlw(LINWTH)
#endif
      END IF
      RETURN

      ENTRY VDSTOS(ATTARR)
      IF (TERMON) THEN
#ifdef Build64
        do 60 i=1,6
          attarr4(i) = attarr(i)
 60     continue
        CALL wx11os(ATTARR4)
#else
        CALL wx11os(ATTARR)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        do 61 i=1,6
          attarr4(i) = attarr(i)
 61     continue
        CALL wmetos(ATTARR4)
#else
        CALL wmetos(ATTARR)
#endif
      END IF
      RETURN

      ENTRY VDWAIT
      IF(TERMON) CALL wx11wt
      IF(FILEON) CALL wmetwt
      RETURN

      ENTRY VDLINA(X,Y)
      IF (TERMON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wx11ln(X4,Y4)
#else
        CALL wx11ln(X,Y)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wmetln(X4,Y4)
#else
        CALL wmetln(X,Y)
#endif
      END IF
      RETURN

      ENTRY VDMOVA(X,Y)
      IF (TERMON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wx11mv(X4,Y4)
#else
        CALL wx11mv(X,Y)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wmetmv(X4,Y4)
#else
        CALL wmetmv(X,Y)
#endif
      END IF
      RETURN

      ENTRY VDNWPG
      IF(TERMON) CALL wx11pg
      IF(FILEON) CALL wmetpg
      RETURN

      ENTRY VDPNTA(X,Y)
      IF (TERMON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wx11pt(X4,Y4)
#else
        CALL wx11pt(X,Y)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        x4 = x
        y4 = y
        CALL wmetpt(X4,Y4)
#else
        CALL wmetpt(X,Y)
#endif
      END IF
      RETURN

      ENTRY VDPOLY(XARRAY,YARRAY,NPTS)
      IF (TERMON) THEN
#ifdef Build64
        if (npts .gt. 1024) stop 'ARRAY BOUNDS IN VDPOLY'
        npts4 = npts
        do 70 i = 1, npts
          xa4(i) = xarray(i)
          ya4(i) = yarray(i)
 70     continue
        CALL wx11py(xa4, ya4, npts4)
#else
        CALL wx11py(XARRAY,YARRAY,NPTS)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        if (npts .gt. 1024) stop 'ARRAY BOUNDS IN VDPOLY'
        npts4 = npts
        do 71 i = 1, npts
          xa4(i) = xarray(i)
          ya4(i) = yarray(i)
 71     continue
        CALL wmetpy(xa4, ya4, npts4)
#else
        CALL wmetpy(XARRAY,YARRAY,NPTS)
#endif
      END IF
      RETURN

      ENTRY VDTERM
      IF(TERMON) CALL wx11tr
      IF(FILEON) CALL wmettr
      RETURN

      ENTRY VDTEXT(LENGTH,CHARS)
      IF (TERMON) THEN
#ifdef Build64
        length4 = length
        if (length .gt. 136) then
          print *, 'Length = ', length
          stop 'ARRAY BOUNDS IN VDTEXT'
        end if
        do 80 i=1, length
          chars4(i) = chars(i)
 80     continue
        CALL wx11tx(LENGTH4,CHARS4)
#else
        CALL wx11tx(LENGTH,CHARS)
#endif
      END IF
      IF(FILEON) THEN
#ifdef Build64
        length4 = length
        if (length .gt. 136) then
          print *, 'Length = ', length
          stop 'ARRAY BOUNDS IN VDTEXT'
        end if
        do 81 i=1, length
          chars4(i) = chars(i)
 81     continue
        CALL wmettx(LENGTH4,CHARS4)
#else
        CALL wmettx(LENGTH,CHARS)
#endif
      END IF
      RETURN

      END
