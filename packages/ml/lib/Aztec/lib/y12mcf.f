      subroutine y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha, aflag,iflag
     *,ifail)
c
c  systens of linear equations are solved by use of sparse matrix tech-
c  nique and by gaussian elimination.
c
      implicit double precision(a-b,g,p,t-y),integer(c,f,h-n,r-s,z)
      double precision a(nn),b(n),pivot(n),aflag(8)
c
c  information which is necessary to begin the elimination is stored.
c
      integer snr(nn),rnr(nn1),ha(iha,11), iflag(10)
      ifail=0
      if(iflag(1).ne.-1)ifail=2
      if(aflag(1).lt.1.0d0)aflag(1)=1.0005 d0
      if(aflag(3).lt.1.0d+5)aflag(3)=1.0d+5
      if(aflag(4).lt.0.0d0)aflag(4)=-aflag(4)
      if(iflag(2).lt.1)ifail=19
      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20
      if(iflag(5).lt.1.or.iflag(5).gt.3)ifail=21
      if(iflag(5).eq.3)ifail=22
      if(ifail.gt.0)go to 1110
      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
c
c  use the information about fill-ins if it is possible.
c
      zz=z
      nr=n*n
      if(iflag(4).ne.2)go to 100
      if(iflag(10).gt.nn)go to 50
      l1=iflag(10)
      l5=l1+1
      if(l5.le.nn)snr(l5)=0
      do 40 i=1,n
      l=n8-i
      l2=ha(l,3)+1
      l3=l2-ha(l,1)
      do 10 j=1,l3
      snr(l5-j)=snr(l2-j)
   10 a(l5-j)=a(l2-j)
      ha(l,3)=l1
      ha(l,1)=l5-l3
      l6=l1-l3
      l5=l5-ha(l,9)
      if(l5.gt.l6)go to 30
      do 20 j=l5,l6
   20 snr(j)=0
   30 continue
   40 l1=l5-1
   50 if(iflag(9).gt.nn1)go to 100
      l2=iflag(9)
      l5=l2+1
      if(l5.le.nn1)rnr(l5)=0
      do 90 i=1,n
      l=n8-i
      l1=ha(l,6)+1
      l4=l1-ha(l,4)
      do 60 j=1,l4
   60 rnr(l5-j)=rnr(l1-j)
      ha(l,4)=l5-l4
      ha(l,6)=l2
      l6=l2-l4
      l5=l5-ha(l,10)
      if(l5.gt.l6)go to 80
      do 70 j=l5,l6
   70 rnr(j)=0
   80 continue
   90 l2=l5-1
  100 r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)
      do 110 i=1,n
      pivot(i)=0.0 d0
      ha(i,2)=ha(i,1)
  110 ha(i,5)=ha(i,4)
      index=ha(n,8)
c
c  start of gaussian elimination.
c
      slut=ha(index,3)-ha(index,2)+1
      do 950 i=1,n7
      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350
      if(iflag(4).ne.2)go to 120
      rrow=ha(i,7)
      rcoll=ha(i,8)
      go to 220
  120 l4=ha(i,8)
      if(iflag(3).eq.1)go to 130
      rrow=l4
      rcoll=rrow
      rpivot=i
      go to 170
  130 r=nr
      v=0.0 d0
      index=iflag(2)
      do 160 kk=1,index
      l1=i-1+kk
      if(l1.gt.n)go to 170
      j=ha(l1,8)
      r7=ha(j,2)
      r8=ha(j,3)
      r9=r8-r7
      t=0.0 d0
      do 140 k=r7,r8
      td=dabs(a(k))
  140 if(t.lt.td)t=td
      t=t/u
      do 160 k=r7,r8
      td=dabs(a(k))
      if(td.lt.t)go to 150
      r6=snr(k)
      r3=r9*(ha(r6,6)-ha(r6,5))
      if(r3.gt.r)go to 150
      if(r3.lt.r)go to 151
      if(v.ge.td)go to 150
  151 v=td
      rrow=j
      rcoll=r6
      r=r3
      rpivot=l1
  150 continue
  160 continue
  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      ha(i,9)=r3
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
c
c  row interchanges.
c
      ha(i,8)=rcoll
  220 if(rrow.eq.i)go to 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t
      do 250 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
      r10=ha(l1,6)
  240 r=r+1
      if(rnr(r).ne.i)go to 240
      rnr(r)=rnr(r10)
  250 rnr(r10)=rrow
      rr3=ha(rrow,2)
      rr4=ha(rrow,3)
      do 270 j=rr3,rr4
      l1=snr(j)
      r=ha(l1,5)-1
  260 r=r+1
      if(rnr(r).ne.rrow)go to 260
  270 rnr(r)=i
      do 280 j=1,3
      r3=ha(rrow,j)
      ha(rrow,j)=ha(i,j)
c
c  column interchanges.
c
  280 ha(i,j)=r3
  290 if(rcoll.eq.i)go to 350
      do 310 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
      r10=ha(l1,3)
  300 r=r+1
      if(snr(r).ne.i)go to 300
      t=a(r10)
      a(r10)=a(r)
      a(r)=t
      snr(r)=snr(r10)
  310 snr(r10)=rcoll
      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)
      do 330 j=c1,cr4
      l1=rnr(j)
      r=ha(l1,2)-1
  320 r=r+1
      if(snr(r).ne.rcoll)go to 320
  330 snr(r)=i
      do 340 j=4,6
      r3=ha(rcoll,j)
      ha(rcoll,j)=ha(i,j)
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  340 ha(i,j)=r3
  350 r9=rr4-rr3
      do 360 rr=rr3,rr4
      if(snr(rr).eq.i)go to 370
  360 continue
      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=dabs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431
      do 430 j=rr3,rr4
      index=snr(j)
  430 pivot(index)=a(j)
  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0)go to 669
      r=rr1
      if(r.gt.rr2)go to 470
      do 460 l=r,rr2
      l1=snr(l)
      td=pivot(l1)
      if(td.eq.0.0d0)go to 450
      pivot(l1)=0.0 d0
      td=a(l)-td*t
      a(l)=td
      td1=dabs(td)
      if(td1.gt.aflag(7))aflag(7)=td1
c
c  too small element is created.remove it from the lists.
c
      if(td1.gt.aflag(2))go to 450
      z=z-1
      a(l)=a(rr1)
      snr(l)=snr(rr1)
      a(rr1)=a(i1)
      snr(rr1)=snr(i1)
      snr(i1)=0
      rr1=rr1+1
      i1=i1+1
      ha(r1,2)=rr1
      ha(r1,1)=i1
      r3=ha(l1,5)
      r2=r3-1
      l4=ha(l1,4)
      l5=rnr(l4)
      l6=rnr(r3)
  440 r2=r2+1
      if(rnr(r2).ne.r1)go to 440
      rnr(r2)=l6
      rnr(r3)=l5
      rnr(l4)=0
      ha(l1,5)=r3+1
      ha(l1,4)=l4+1
  450 continue
  460 continue
  470 continue
      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0d0)go to 740
      tol3=-tol2*t
      tol1=dabs(tol3)
      if(tol1.lt.aflag(2))go to 740
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2
      if(iflag(4).ne.1)go to 480
      if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
      if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
  480 if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
c
c  collection in row ordered list.
c
  500 r10=nn-lfr
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1
      do 520 jj=1,n
      l1=ha(jj,3)
      if(l1.lt.ha(jj,1))go to 510
      ha(jj,3)=snr(l1)
      snr(l1)=-jj
  510 continue
  520 continue
      l3=0
      l4=1
      do 550 jj=1,r4
      if(snr(jj).eq.0)go to 540
      l3=l3+1
      if(snr(jj).gt.0)go to 530
      l5=-snr(jj)
      snr(jj)=ha(l5,3)
      ha(l5,3)=l3
      l6=l4+ha(l5,2)-ha(l5,1)
      ha(l5,2)=l6
      ha(l5,1)=l4
      l4=l3+1
  530 a(l3)=a(jj)
      snr(l3)=snr(jj)
  540 continue
  550 continue
      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j
      if(r10.ge.r4)go to 560
      ifail=5
c
c fill-in takes place in the row ordered list.
c
      go to 1110
  560 r8=lfr-1
      rr2=r4+lfr
      if(r8.le.0)go to 579
      l3=i1-1
      do 570 ll=1,r8
      l4=r4+ll
      l5=l3+ll
      a(l4)=a(l5)
      snr(l4)=snr(l5)
  570 snr(l5)=0
  579 rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590
  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=dabs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z) iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
c
c  collection in column ordered list.
c
  630 r10=nn1-lfc
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1
      do 640 jj=i,n
      l1=ha(jj,6)
      ha(jj,6)=rnr(l1)
  640 rnr(l1)=-jj
      l3=0
      l4=1
      do 670 jj=1,r5
      if(rnr(jj).eq.0)go to 660
      l3=l3+1
      if(rnr(jj).gt.0)go to 650
      l5=-rnr(jj)
      rnr(jj)=ha(l5,6)
      ha(l5,6)=l3
      l6=l4+ha(l5,5)-ha(l5,4)
      ha(l5,5)=l6
      ha(l5,4)=l4
      l4=l3+1
  650 rnr(l3)=rnr(jj)
  660 continue
  670 continue
      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)
      if(r10.ge.r5)go to 680
      ifail=6
c
c fill-in takes place in the column ordered list.
c
      go to 1110
  680 r8=lfc-1
      cr2=r5+lfc
      if(r8.le.0)go to 699
      l3=c2-1
      do 690 l=1,r8
      l4=r5+l
      l5=l3+l
      rnr(l4)=rnr(l5)
  690 rnr(l5)=0
  699 cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710
  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
      go to 1110
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue
      if(r9.le.0)go to 882
      do 881 j=rr3,rr4
      index=snr(j)
  881 pivot(index)=0.0 d0
  882 continue
      cr3=ha(i,4)
      do 890 j=cr3,cr4
  890 rnr(j)=0
      if(r9.le.0)go to 930
      l2=ha(i,2)-1
      do 920 ll=1,r9
      r=snr(l2+ll)
      r1=ha(r,5)
      r2=ha(r,6)
      if(r2.gt.r1)go to 900
      ifail=8
      go to 1110
  900 ha(r,5)=r1+1
      r3=r1-1
  910 r3=r3+1
      if(rnr(r3).ne.i)go to 910
      rnr(r3)=rnr(r1)
  920 rnr(r1)=i
  930 aflag(5)=aflag(7)/aflag(6)
      if(aflag(5).lt.aflag(3))go to 940
      ifail=4
      go to 1110
  940 continue
c
c  preparation to begin the back substitution.
c
  950 continue
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0 d0
      td=dabs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td
      if(td.gt.grmin)go to 960
      ifail=3
      go to 1110
  960 if(iflag(4).ne.1)go to 1060
      iflag(10)=ha(n,9)
      iflag(9)=ha(n,10)
      do 990 i=1,n7
      r1=n-i
      iflag(10)=iflag(10)+ha(r1,9)
      iflag(9)=iflag(9)+ha(r1,10)
      if(iflag(3).eq.0)go to 980
      do 970 j=9,10
      r2=ha(r1,j-2)
      r6=ha(r2,j)
      ha(r2,j)=ha(r1,j)
  970 ha(r1,j)=r6
  980 continue
  990 continue
1060  continue
      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end
