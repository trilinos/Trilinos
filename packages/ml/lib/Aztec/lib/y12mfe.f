      subroutine y12mfe(n,a,snr,nn,rnr,nn1,a1,sn,nz,ha,iha,b,b1,x,y,
     1 aflag,iflag,ifail)
      implicit real(a-b,d,g,p,t-y),integer(c,f,h-n,r-s,z)
      double precision ext,er,er1,er2,err,e
      real a(nn),b(n),b1(n),x(n),y(n),a1(nz),aflag(11)
      integer snr(nn),rnr(nn1),ha(iha,13),sn(nz),iflag(12)
c
c  store the non-zero elements,their column numbers,information about
c  row starts,information about row ends and the right-hand side.
c
      ifail=0
      nres=0
      dres=0.0
      state=iflag(5)
      kit=1
      it=iflag(11)
      if(state.eq.1)ifail=10
      if(it.lt.2)ifail=23
      if(ifail.ne.0)go to 160
      do 10 i=1,n
   10 b1(i)=b(i)
      if(state.eq.3)go to 70
      call y12mbe(n,nz,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 160
      do 20 i=1,nz
      sn(i)=snr(i)
   20 a1(i)=a(i)
      do 30 i=1,n
      ha(i,12)=ha(i,1)
   30 ha(i,13)=ha(i,3)
      if(aflag(2).ge.0.0)go to 60
      gt1=aflag(6)
      do 50 i=1,n
      l1=ha(i,1)
      l2=ha(i,3)
      gt2=0.0
      do 40 j=l1,l2
      d=abs(a(j))
   40 if(gt2.lt.d)gt2=d
   50 if(gt2.lt.gt1)gt1=gt2
      aflag(2)=-gt1*aflag(2)
c
c  find the first solution.
c
   60 call y12mce(n,nz,a,snr,nn,rnr,nn1,y,b,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 160
   70 call y12mde(n,a,nn,b,y,snr,ha,iha,iflag,ifail)
      if(ifail.ne.0)go to 160
c
c  prepare the data in order to begin the iterations.
c
      dd=0.0
      do 80 i=1,n
      x(i)=b(i)
      xx=abs(b(i))
   80 if(dd.lt.xx)dd=xx
      xm=dd
      if(dd.eq.0.)go to 160
c
c  begin to iterate.
c
   90 d=dd
      dres=0.
      do 110 i=1,n
      er=b1(i)
      l1=ha(i,12)
      l2=ha(i,13)
      do 100 j=l1,l2
      er1=a1(j)
      l3=sn(j)
      er2=x(l3)
  100 er=er-er1*er2
c
c  store residuals rounded to single precision.
c
      b(i)=er
      xx=dabs(er)
  110 if(dres.lt.xx)dres=xx
      if(dres.eq.0.)go to 160
      if(nres.eq.1) go to 150
      if(dres.gt.1.0e+4*xm)go to 150
      kit=kit+1
      iflag(5)=3
      call y12mde(n,a,nn,b,y,snr,ha,iha,iflag,ifail)
      if(ifail.ne.0)go to 160
c
c  compute the uniform norm of the current solution vector.
c
      dd=0.
      do 120 i=1,n
      xx=abs(b(i))
  120 if(dd.lt.xx)dd=xx
      if(dd.eq.0.0)go to 160
c
c  check the convergence criterion.
c
      if(dd.gt.d.and.kit.gt.2)go to 160
c
c  calculate an improved solution.
c
      dcor=dd
      xm=0.0
      do 130 i=1,n
      x(i)=x(i)+b(i)
      xx=abs(x(i))
  130 if(xx.gt.xm) xm=xx
c
c  check the stopping criteria.
c
      if(10.0+dd/xm.eq.10.0) go to 140
      if(kit.lt.it) go to 90
c
c  end of the iterations.
c
  140 nres=1
      go to 90
  150 dd=abs(dd)
  160 iflag(5)=state
      iflag(12)=kit
      aflag(9)=dd
      aflag(10)=dres
      aflag(11)=xm
      return
      end
