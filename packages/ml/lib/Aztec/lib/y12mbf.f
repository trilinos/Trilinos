      subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag,
     1 iflag,ifail)
c
c
c  the non-zero elements of a sparse matrix a are prepared  in order to
c  solve the system ax=b by use of sparse matrix technique/
c
c
      implicit double precision(a-b,g,p,t-y),integer(c,f,h-n,r-s,z)
      double precision a(nn), aflag(8)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      mode=iflag(4)
      ifail=0
      if(n.lt.2)ifail=12
      if(z.le.0)ifail=13
      if(nn.lt.2*z)ifail=5
      if(nn1.lt.z)ifail=6
      if(ifail.eq.0.and.n.gt.z)ifail=14
      if(iha.lt.n)ifail=15
      if(mode.lt.0)ifail=16
      if(mode.gt.2)ifail=16
      if(ifail.ne.0) go to 22
      gt1=0.0d0
      do 10 i=1,n
      ha(i,2)=0
      ha(i,3)=0
   10 ha(i,6)=0
c
c  find the number of the non-zero elements in each row and column;move
c  the non-zero elements in the end of the arrays a and snr;find the
c  largest non-zero element in a(in absolute value).
c
      do 20 i=1,z
      t=dabs(a(i))
      l3=rnr(i)
      l4=snr(i)
      if(l4.gt.n.or.l4.lt.1)ifail=24
      if(l3.gt.n.or.l3.lt.1)ifail=25
      ha(l3,3)=ha(l3,3)+1
      ha(l4,6)=ha(l4,6)+1
      if(t.gt.gt1)gt1=t
      a(z+i)=a(i)
   20 snr(z+i)=snr(i)
      if(ifail.gt.0)go to 22
c
c  store the information of the row starts(in ha(i,1))and of the column
c  starts(in ha(i,4)).
c
      l1=1
      l2=1
      do 40 i=1,n
      l3=ha(i,3)
      l4=ha(i,6)
      if(l3.gt.0)go to 21
      ifail=17
      go to 22
   21 if(l4.gt.0)go to 23
      ifail=18
      go to 22
   23 if(mode.eq.2)go to 30
      ha(i,9)=l3
      ha(i,10)=l4
      ha(i,11)=0
      ha(l3,2)=ha(l3,2)+1
      ha(i,5)=l3
   30 ha(i,1)=l1
      ha(i,4)=l2
      l1=l1+l3
      l2=l2+l4
      ha(i,3)=0
   40 ha(i,6)=0
c
c  store the non-zero elements of matrix a(ordered in rows) in the
c  first z locations of the array a.do the same for their column numbers
c
      do 50 i=1,z
      l1=z+i
      l3=rnr(i)
      l2=ha(l3,1)+ha(l3,3)
      a(l2)=a(l1)
      snr(l2)=snr(l1)
   50 ha(l3,3)=ha(l3,3)+1
c
c  store the row numbers of the non-zero elements ordered by columns in
c  the first z locations of the array rnr. store information about row
c  ends(in ha(i,3)).
c
      l4=1
      do 70 i=1,n
      if(mode.eq.2)go to 60
      if(ha(i,2).eq.0)go to 60
      ha(i,11)=l4
      l4=l4+ha(i,2)
      ha(i,2)=ha(i,11)
   60 ha(i,3)=ha(i,1)+ha(i,3)-1
      l1=ha(i,1)
      l2=ha(i,3)
      do 70 j=l1,l2
      l3=snr(j)
      r=ha(l3,6)
      index=ha(l3,4)+r
      rnr(index)=i
      if(r.eq.0)go to 70
      if(j.eq.l1)go to 70
      if(rnr(index-1).ne.i)go to 70
      ifail=11
      go to 22
   70 ha(l3,6)=r+1
      do 90 i=1,n
      if(mode.eq.2)go to 80
      l3=ha(i,5)
      l5=ha(l3,2)
      ha(l5,8)=i
      ha(i,7)=l5
      ha(l3,2)=ha(l3,2)+1
   80 continue
   90 ha(i,6)=ha(i,4)+ha(i,6)-1
      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1
22    return
      end
