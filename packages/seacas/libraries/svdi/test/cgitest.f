C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c      external cgit07
c      external cgimet
      character*10 met1,met2
      integer mf1id, tk1id
      real x(20), y(20)
      integer clr(3)
      integer colors(24)
      character*10 text
      integer lnclr(3)
      integer chhit(10)

      include 'cgidef.f'

c set colors up like vdi:
C           Index  Color  RGB Values
C             0    black   0.,0.,0.
C             1    red     1.,0.,0.
C             2    green   0.,1.,0.
C             3    yellow  1.,1.,0.
C             4    blue    0.,0.,1.
C             5    magenta 1.,0.,1.
C             6    cyan    0.,1.,1.
C             7    white   1.,1.,1.

      data colors / 0, 0, 0,  255, 0, 0, 0, 255, 0,
     *               255, 255, 0, 0, 0, 255,  255, 0, 255,
     *               0, 255, 255,  255, 255, 255 /

      met1 = 'file55'
      met2 = 'meta2'

c      call xcact(cgimet,mf1id)
c      call xcact(cgimet,mf2id)
c      call xcact(cgit07,tk1id)

c      call xcooff(tk1id)
c      call xcooff(mf2id)
c      call cesc(-28372,1,met1)
c      call xcooff(mf1id)
c      call xcoon(mf2id)
c      call cesc(-28372,1,met2)
c      call xcoon(tk1id)
c      call xcoon(mf1id)

c initialize cgi
c      call cesc(-28372,1,met1)
      call ci(1)
      call cvdcx(0., 0., 32767.,32767.)
c      call cv(.5,0.,2.5,1.)
C VDC clip COFF, CON
      call ccl(CON)
C viewsurface clip CDCOFF, CDCREC, CVPORT
      call cdscl(CDCREC)

      do 10 i=1,2
c draw some primitives in default color
c  ..polylines
      x(1) = 100.
      x(2) = 16000.
      y(1) = 16500.
      y(2) = 16500.
      call clnt(1)
      call cpl(2,x,y)
      y(1) = 20500.
      y(2) = 20500.
      call clnt(2)
      call cpl(2,x,y)
      y(1) = 24500.
      y(2) = 24500.
      call clnt(3)
      call cpl(2,x,y)
      y(1) = 28500.
      y(2) = 28500.
      call clnt(4)
      call cpl(2,x,y)

c  ..polygons
      call cis(CHOLLO)
      x(1) = 20000.
      y(1) = 25000.
      x(2) = 29000.
      y(2) = 25000.
      x(3) = 24500.
      y(3) = 32000
      x(4) = 20000.
      y(4) = 25000.
      call cpg(4,x,y)
      x(1) = 20000.
      y(1) = 23500.
      x(2) = 29000.
      y(2) = 23500.
      x(3) = 24500.
      y(3) = 16500.
      x(4) = 20000.
      y(4) = 23500.
      call cis(CSOLID)
      call cpg(4,x,y)

c ...markers
      x(1) = 100.
      y(1) = 100.
      x(2) = 2100.
      y(2) = 2100.
      x(3) = 4100.
      y(3) = 4100.
      x(4) = 6100.
      y(4) = 6100
      x(5) = 8100.
      y(5) = 8100.
      x(6) = 10100.
      y(6) = 10100.
      x(7) = 12100.
      y(7) = 12100.
      x(8) = 14100.
      y(8) = 14100.
      x(9) = 16100.
      y(9) = 16100.
      call cpm(9,x,y)

c ...text
c ...do tektronix first, then metafile
c ...do this because cgtxx is an inquiry
c      do 20 j=1,2
c      if(j.eq.1)then
c        call xcooff( mf1id )
c        call xcooff( mf2id )
c        call xcoon( tk1id )
c        call xcsol( tk1id )
c        call cili( CLOCAT, 1)
c      else
c        call xcooff( tk1id )
c        call xcoon( mf1id )
c        call xcoon( mf2id )
c        call xcsol( mf1id )
        call cili( CLOCAT, 1)
c     endif
      call cqchh(" ",0,1,1,istat,ntot,nlist,chhit)
      text = "TEXT"
      call cchh(327.)
      call ctx( 17000.,100.,1,text(1:4))
      call cgtxx(17000.,100.,text(1:4),ivstat,ivconc,xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      call ctx( xconc,yconc,1,text(1:4))
      call cgtxx(xconc,yconc,text(1:4),ivstat,ivconc, xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      call cchh(3270.)
      call ctx( xconc,yconc,1,text(1:4))

      call cchh(3270.)
      text = "A"
      call ctx( 16500.,10000.,1,text(1:1))
      call cgtxx(16500.,10000.,text(1:1),ivstat,ivconc, xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      text = "B"
      call ctx( xconc,yconc,1,text(1:1))
      call cgtxx(xconc,yconc,text(1:1),ivstat,ivconc, xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      text = "C"
      call ctx( xconc,yconc,1,text(1:1))
      call cgtxx(xconc,yconc,text(1:1),ivstat,ivconc, xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      text = "ABC"
      call ctx( xconc, yconc,1,text(1:3))
      call cgtxx(xconc,yconc,text(1:3),ivstat,ivconc, xconc,yconc,
     *  p1,q1,p2,q2,p3,q3,p4,q4)
      x(1) = p1
      y(1) = q1
      x(2) = p2
      y(2) = q2
      x(3) = p3
      y(3) = q3
      x(4) = p4
      y(4) = q4
      x(5) = p1
      y(5) = q1
      call clnt(1)
      call cpl(5,x,y)

c use request locator as a pause
      call crqlc(1,1,istat,irstat,mvalid,itrig,xx,yy)
      call cpds(0)
20    continue

c turn everything on
c      call xcoon( tk1id )
c      call xcoon( mf1id )
c      call xcoon( mf2id )

c change colors
      if(i.eq.1) then

c does the device support direct color
c ...do it this weird way cause i'm lazy - if one is only indexed
c     color, then treat them both as indexed color
c        call xcsol( tk1id )
        call cqc( idum, idum, idum, idum, icmode, idum, idum, idum )
        if( icmode .eq. CCLRID ) then
c          call xcsol( mf1id )
          call cqc( idum, idum, idum, idum, icmode, idum, idum, idum )
        endif

c set direct color mode if supported, otherwise use indexed
c only set the color table for indexed color
        if( icmode .eq. CCLRID ) then
          call ccsm( CDRECT )
        else
          call ccsm( CINDEX )
          call cct(0,8,colors)
        endif

c  ...make polylines red
        if( icmode .eq. CCLRID ) then
          clr(1) = 255
          clr(2) = 0
          clr(3) = 0
        else
          clr(1) = 1
        endif
        call clnc(clr)
        call cqlna(istat,lnbi,lntyp,lwmod,lnwid,csmod,lnclr,lcmod)

c  ...make polygons green
        if( icmode .eq. CCLRID ) then
          clr(1) = 0
          clr(2) = 255
          clr(3) = 0
        else
          clr(1) = 2
        endif
        call cflc(clr)

c  ...make polymarkers yellow
        if( icmode .eq. CCLRID ) then
          clr(1) = 255
          clr(2) = 255
          clr(3) = 0
        else
          clr(1) = 3
        endif
        call cmkc(clr)

c  ...make text blue
        if( icmode .eq. CCLRID ) then
          clr(1) = 0
          clr(2) = 0
          clr(3) = 255
        else
          clr(1) = 4
        endif
        call ctxc(clr)
      endif

10    continue
      call ct
c      call xcdact( mf1id )
c      call xcdact( mf2id )
c      call xcdact( tk1id )
      stop
      end
