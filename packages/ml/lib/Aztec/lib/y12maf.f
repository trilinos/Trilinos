      subroutine y12maf(n, z, a, snr, nn, rnr, nn1, pivot, ha,
     1iha,aflag,iflag,b,ifail)
      implicit double precision (a-b,g,p,t-y), integer (c,f,h-n,r-s,z)
      double precision a(nn), pivot(n), aflag(8),b(n)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      aflag(1)=16.0d0
      aflag(2)=1.d-12
      aflag(3)=1.d+16
      aflag(4)=1.d-12
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1
      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 1
      call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     1 ifail)
      if(ifail.ne.0)go to 1
      call y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
1     return
      end
