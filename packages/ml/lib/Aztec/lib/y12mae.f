      subroutine y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha,
     1iha,aflag,iflag,b,ifail)
      implicit real (a-b,g,p,t-y), integer (c,f,h-n,r-s,z)
      real a(nn), pivot(n), aflag(8),b(n)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      aflag(1)=16.
      aflag(2)=1.e-12
      aflag(3)=1.e+16
      aflag(4)=1.e-12
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1
      call y12mbe(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 1
      call y12mce(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     1 ifail)
      if(ifail.ne.0)go to 1
      call y12mde(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
1     return
      end
