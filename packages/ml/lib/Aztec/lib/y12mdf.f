      subroutine y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
      implicit double precision(a-b,g,p,t-y),integer (c,f,h-n,r-s,z)
      double precision a(nn), pivot(n), b(n)
      integer snr(nn), ha(iha,11), iflag(10)
      ifail=0
      if(iflag(1).eq.-2)go to 1000
      ifail=1
      go to 1110
1000  mode=iflag(4)
      ipiv=iflag(3)
      n8=n+1
      n7=n-1
      state=iflag(5)
c
c  solve the system with lower triangular matrix  l  (if the
c  lu-factorization is available).
c
      if(state.ne.3)go to 1051
      if(ipiv.eq.0)go to 1020
      do 1010 i=1,n7
      l1=ha(i,7)
      t=b(l1)
      b(l1)=b(i)
      b(i)=t
1010  continue
1020  continue
      do 1050 i=1,n
      rr1=ha(i,1)
      rr2=ha(i,2)-1
      if(rr1.gt.rr2)go to 1040
      do 1030 j=rr1,rr2
      l1=snr(j)
 1030 b(i)=b(i)-a(j)*b(l1)
 1040 continue
 1050 continue
c
c  solve the system with upper triagular matrix.
c
 1051 continue
      do 1090 i=1,n
      r1=n8-i
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      if(rr2.lt.rr1)   go to 1080
      do 1070 j=rr1,rr2
      r2=snr(j)
 1070 b(r1)=b(r1)-a(j)*b(r2)
 1080 continue
 1090 b(r1)=b(r1)/pivot(r1)
c
c if interchanges were used during the  elimination then a reordering in
c lution vector is made.
c
      if(ipiv.eq.0)go to 1110
      do 1100 i=1,n7
      r1=n-i
      r2=ha(r1,8)
      t=b(r2)
      b(r2)=b(r1)
 1100 b(r1)=t
 1110 return
      end
