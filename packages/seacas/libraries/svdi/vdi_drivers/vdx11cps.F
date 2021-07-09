C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C   cps driver is the file device
C   x11 driver is the terminal device
      SUBROUTINE VDESCP(ESCPCD,N,ARGS)
      INTEGER ESCPCD,SUPPRT,TRMSUP,FILSUP
      REAL ARGS(*)
      LOGICAL TERMON,FILEON,TERMUP,FILEUP

      real locx, locy
      real clrary(*)
      integer clrmod
      integer indexa(*)
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
      integer*4 trmsup4, supprt4, esc4
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
#else
#endif
      DATA TERMON,FILEON,TERMUP,FILEUP /2*.TRUE.,2*.TRUE./

      IF ( ESCPCD.GE.10000 .AND. ESCPCD.LE.10003 ) THEN

C Special Escape Codes for Dual Device Control:
         KSTAT = ESCPCD - 10000
         TERMON = (KSTAT.EQ.1 .OR. KSTAT.EQ.3) .AND. TERMUP
         FILEON = (KSTAT.EQ.2 .OR. KSTAT.EQ.3) .AND. FILEUP
      ELSE
         IF (TERMON) THEN
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
         IF (FILEON) CALL wcpses(ESCPCD,N,ARGS)
      END IF
      RETURN

C------------------------------------------------------------------------
      ENTRY VDIQES(ESCPCD,SUPPRT)

C Special Escape Codes for Dual Device Control:
      IF ( ESCPCD.EQ.10000 ) THEN
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
         CALL wcpsie(10000,FILSUP)
         SUPPRT = 1 - FILSUP
      ELSE IF ( ESCPCD.EQ.10003 ) THEN
#ifdef Build64
         esc4 = 10000
         CALL wx11ie(esc4,TRMSUP4)
         TRMSUP = TRMSUP4
#else
         CALL wx11ie(10000,TRMSUP)
#endif
         CALL wcpsie(10000,FILSUP)
         SUPPRT = (1 - TRMSUP) * (1 - FILSUP)
      ELSE
         IF (FILEON) CALL wcpsie(ESCPCD,SUPPRT)
         IF (TERMON) THEN
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

C-----------------------------------------------------------------------
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
      CALL wcpsie( 10000,ISTAT )
      FILEUP = ISTAT .EQ. 0
      FILEON = FILEON .AND. FILEUP

C Stop program if both devices are null.
      IF (.NOT.(TERMUP.OR.FILEUP)) STOP 'No graphics device assigned.'

      IF (TERMON) THEN
#ifdef Build64
         ASPECT4 = ASPECT
         JUSTIF4 = JUSTIF
         CALL wx11nt(ASPECT4,JUSTIF4)
#else
         CALL wx11nt(ASPECT,JUSTIF)
#endif
      END IF
      IF (FILEON) CALL wcpsnt(ASPECT,JUSTIF)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDAABU(BTNNUM)
      IF (FILEON) CALL wcpsbu(BTNNUM)
      IF (TERMON) THEN
#ifdef Build64
C ... Not implemented in X11 driver
         CALL wx11bu(BTNNUM)
#else
         CALL wx11bu(BTNNUM)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDABGL(BTNNUM,X,Y)
      IF (FILEON) CALL wcpsbl(BTNNUM,X,Y)
      IF (TERMON) THEN
#ifdef Build64
C ... Not implemented in X11 driver
         CALL wx11bl(BTNNUM,X,Y)
#else
         CALL wx11bl(BTNNUM,X,Y)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDAKGL(CHAR,X,Y)
      IF (FILEON) CALL wcpskl(CHAR,X,Y)
      IF (TERMON) THEN
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

C-----------------------------------------------------------------------
      ENTRY VDALOC(X,Y)
      IF (FILEON) CALL wcpslo(X,Y)
      IF (TERMON) THEN
#ifdef Build64
         CALL wx11lo(X4,Y4)
         x = x4
         y = y4
#else
         CALL wx11lo(X,Y)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDBELL
      IF (TERMON) THEN
         CALL wx11be
      END IF
      IF (FILEON) CALL wcpsbe
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDBUFL
      IF (TERMON) THEN
         CALL wx11fl
      END IF
      IF (FILEON) CALL wcpsfl
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDFRAM(ITYPE)
      IF (TERMON) THEN
C ... Doesn't do anything in X11 driver
         CALL wx11fr(ITYPE)
      END IF
      IF (FILEON) CALL wcpsfr(ITYPE)
      RETURN

C .... disable below
C-----------------------------------------------------------------------
      ENTRY VDIQCO(NUM,INDEXA,CLRARY,CLRMOD)
      IF (FILEON) CALL wcpsic(NUM,INDEXA,CLRARY,CLRMOD)
      IF (TERMON) THEN
#ifdef Build64
         num4 = num
         clrmod4 = clrmod
         do 20 i=1, num
            index4a(i) = indexa(i)
 20      continue
         CALL wx11ic(NUM4,INDEX4a,CLRARY4,CLRMOD4)
         do 30 i=1, 3*num
            clrary(i) = clrary4(i)
 30      continue
#else
         CALL wx11ic(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDIQCP(X,Y)
      IF (FILEON) CALL wcpscp(X,Y)
      IF (TERMON) THEN
#ifdef Build64
         CALL wx11cp(X4,Y4)
         x = x4
         y = y4
#else
         CALL wx11cp(X,Y)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDIQDC(INDEX,VALUE)
      IF (FILEON) CALL wcpsdc(INDEX,VALUE)
      IF (TERMON) THEN
#ifdef Build64
         index4 = index
         CALL wx11dc(INDEX4,VALUE4)
         value = value4
#else
         CALL wx11dc(INDEX,VALUE)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDIQND(XNDC,YNDC)
      IF (FILEON) CALL wcpsid(XNDC,YNDC)
      IF (TERMON) THEN
#ifdef Build64
        CALL wx11id(XNDC4,YNDC4)
        xndc = xndc4
        yndc = yndc4
#else
        CALL wx11id(XNDC,YNDC)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDIQOS(ATTARR)
      IF (FILEON) CALL wcpsio(ATTARR)
      IF (TERMON) THEN
#ifdef Build64
         CALL wx11io(ATTARR4)
         do 40 i=1,7
            attarr(i) = attarr4(i)
 40      continue
#else
         CALL wx11io(ATTARR)
#endif
      END IF
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDLOGE(ERRNUM,ERRSEV)
      IF (TERMON) THEN
C ... Doesn't do anything in X11 driver
         CALL wx11le(ERRNUM,ERRSEV)
      END IF
      IF (FILEON) CALL wcpsle(ERRNUM,ERRSEV)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTBC(COLOR)
      IF (TERMON) THEN
#ifdef Build64
         color4 = color
         CALL wx11bc(COLOR4)
#else
         CALL wx11bc(COLOR)
#endif
      END IF
      IF (FILEON) CALL wcpsbc(COLOR)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTCO(NUM,INDEXA,CLRARY,CLRMOD)

      IF (TERMON) THEN
#ifdef Build64
         num4 = num
         clrmod4 = clrmod
         j = 1
         do 50 i=1, num
            index4a(i) = indexa(i)
 50      continue
         do 55 i=1, num*3
            clrary4(i) = clrary(i)
 55       continue
         CALL wx11co(NUM4,INDEX4A,CLRARY4,CLRMOD4)
#else
         CALL wx11co(NUM,INDEXA,CLRARY,CLRMOD)
#endif
      END IF
      IF (FILEON) CALL wcpsco(INDEXA,CLRARY,CLRMOD)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTCS(YSIZE)
      IF (TERMON) THEN
#ifdef Build64
         Y4 = YSIZE
         CALL wx11cs(Y4)
#else
         CALL wx11cs(YSIZE)
#endif
      END IF
      IF (FILEON) CALL wcpscs(YSIZE)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTFC(COLOR)
      IF (TERMON) THEN
#ifdef Build64
         color4 = color
         CALL wx11fc(COLOR4)
#else
         CALL wx11fc(COLOR)
#endif
      END IF
      IF (FILEON) CALL wcpsfc(COLOR)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTIN(INTEN)
      IF (TERMON) THEN
#ifdef Build64
         inten4 = inten
         CALL wx11in(INTEN4)
#else
         CALL wx11in(INTEN)
#endif
      END IF
      IF (FILEON) CALL wcpsin(INTEN)
      RETURN

C-----------------------------------------------------------------------
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
      IF (FILEON) CALL wcpsla(LOCX,LOCY)
      RETURN

C....enabled below
C-----------------------------------------------------------------------
      ENTRY VDSTLS(LINSTY)
      IF (TERMON) THEN
#ifdef Build64
         linsty4 = linsty
         CALL wx11ls(LINSTY4)
#else
         CALL wx11ls(LINSTY)
#endif
      END IF
      IF (FILEON) CALL wcpsls(LINSTY)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTLW(LINWTH)
      IF (TERMON) THEN
#ifdef Build64
         linwth4 = linwth
         CALL wx11lw(LINWTH4)
#else
         CALL wx11lw(LINWTH)
#endif
      END IF
      IF (FILEON) CALL wcpslw(LINWTH)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDSTOS(ATTARR)
      IF (TERMON) THEN
#ifdef Build64
         do 60 i=1,6
            attarr4(i) = attarr(i)
 60      continue
         CALL wx11os(ATTARR4)
#else
         CALL wx11os(ATTARR)
#endif
      END IF
      IF (FILEON) CALL wcpsos(ATTARR)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDWAIT
      IF (TERMON) THEN
         CALL wx11wt
      END IF
      IF (FILEON) CALL wcpswt
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDLINA(XX,YY)
      IF (TERMON) THEN
#ifdef Build64
         x4 = xx
         y4 = yy
         CALL wx11ln(X4,Y4)
#else
         CALL wx11ln(XX,YY)
#endif
      END IF
      IF (FILEON) CALL wcpsln(XX,YY)
      RETURN

C-----------------------------------------------------------------------
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
      IF (FILEON) CALL wcpsmv(X,Y)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDNWPG
      IF (TERMON) THEN
         CALL wx11pg
      END IF
      IF (FILEON) CALL wcpspg
      RETURN

C-----------------------------------------------------------------------
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
      IF (FILEON) CALL wcpspt(X,Y)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDPOLY(XARRAY,YARRAY,NPTS)
      IF (TERMON) THEN
#ifdef Build64
         if (npts .gt. 1024) stop 'ARRAY BOUNDS IN VDPOLY'
         npts4 = npts
         do 70 i = 1, npts
            xa4(i) = xarray(i)
            ya4(i) = yarray(i)
 70      continue
         CALL wx11py(xa4, ya4, npts4)
#else
         CALL wx11py(XARRAY,YARRAY,NPTS)
#endif
      END IF
      IF (FILEON) CALL wcpspy(XARRAY,YARRAY,NPTS)
      RETURN

C-----------------------------------------------------------------------
      ENTRY VDTERM
      IF (TERMON) THEN
         CALL wx11tr
      END IF
      IF (FILEON) CALL wcpstr
      RETURN

C-----------------------------------------------------------------------
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
 80      continue
         CALL wx11tx(LENGTH4,CHARS4)
#else
         CALL wx11tx(LENGTH,CHARS)
#endif
      END IF
      IF (FILEON) CALL wcpstx(LENGTH,CHARS)
      RETURN

      END

      SUBROUTINE VBERRH()
      RETURN
      END
