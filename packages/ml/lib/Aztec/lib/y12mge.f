      subroutine y12mge(n,nn,a,snr,w,pivot,anorm,rcond,iha,ha
     *                 ,iflag,ifail)
c
c
c   purpose.
c   --------
c
c   this subroutine computes the number called    rcond    by
c   dongarra et al.(1979). this number is the reciprocal of the
c   condition number of matrix   a .   the subroutine can be
c   called immediately after the call of   y12mce.   the parameters
c   (except   rcond   and   anorm ) have the same meaning as the
c   corresponding parameters in the subroutines of package   y12m  (the
c   subroutine can be  call only if the   lu   decomposition of matrix
c   a   computed by   y12mce   is not destroyed). subroutine  y12mhe
c   should be called before the call of subroutine   y12mce (this
c   subroutine calculates the   one-norm   of matrix   a   and
c   stores  it  in   anorm.   on successful exit   rcond   will
c   contain an approximation to the reciprocal of the condition
c   number of matrix   a.   more details, concerning the use
c   of similar subroutines for  dense matrices, can be found
c   in   j.dongarra, j.r.bunch, c.b.moler and g.w.stewart (1979):
c        "linpack - user's guide", siam, philadelphia.
c
c
c
c
c  declaration of the global variables and arrays.
c
c
      integer n, nn, iha, iflag, ifail
      integer   snr, ha
      real      a, w,    pivot, rcond, anorm
      dimension a(nn),snr(nn),ha(iha,3),pivot(n),w(n),iflag(5)
c
c
c  declaration of the internal variables.
c
c
      real  aa, ynorm, znorm, t
      integer  l1, l2, l3, l, n7, n8
c
c
c   check whether the entry is correct or not.
c
c
      if(ifail.ne.0) go to 180
      if(iflag(5).ne.1) go to 10
      ifail=26
      go to 180
c
c
c   no error detected.  the computations will be continued.
c
c
   10 n8=n+1
      n7=n-1
c
c
c   solve a system of the form    u1*w=e   where   u1   is the
c   transpose of matrix   u   in the   lu-factorization of matrix
c   a   and    e    is a vector whose components are equal to   +1
c   or   -1.
c
c
      w(1)=1.0/pivot(1)
      do   20   i=2,n
   20 w(i)=0.0
      do   50   i=2,n
      l1=ha(i,2)
      l2=ha(i,3)
      if(l1.gt.l2) go to 40
      t=w(i-1)
      do   30   j=l1,l2
      l=snr(j)
   30 w(l)=w(l)+t*a(j)
   40 if(w(i).gt.0.0) w(i)=w(i)+1.0
      if(w(i).le.0.0) w(i)=w(i)-1.0
   50 w(i)=-w(i)/pivot(i)
c
c
c   solve a system of the form   l1*y=w   where   l1   is the
c   transpose of matrix   l   in the   lu-factorization of
c   matrix   a .   the components of vector   y   are stored
c   array   w  (thus, the contents of array   w   are overwritten
c   by the components of vector   y ).
c
c
      do   80   i=1,n7
      l=n-i
      l1=ha(l,1)
      l2=ha(l,2)-1
      if(l1.gt.l2) go to 70
      t=w(l+1)
      do   60   j=l1,l2
      l3=snr(j)
   60 w(l3)=w(l3)-t*a(j)
   70 continue
   80 continue
c
c
c   calculate the one-norm of vector   y .
c
c
      ynorm=0.0
      do   90   i=1,n
   90 ynorm=ynorm+abs(w(i))
c
c
c   compute the solution of    (lu)z=y .  this means that
c   two systems with triangular matrices are solved using the
c   same ideas as above. the components of the calculated solution
c   are stored in array   w .
c
c
      do 130 i=1,n
      l1=ha(i,1)
      l2=ha(i,2)-1
      if(l1.gt.l2) go to 120
      do 110 j=l1,l2
      l=snr(j)
  110 w(i)=w(i)-a(j)*w(l)
  120 continue
  130 continue
      do 160 i=1,n
      l3=n8-i
      l1=ha(l3,2)
      l2=ha(l3,3)
      if(l1.gt.l2) go to 150
      do 140 j=l1,l2
      l=snr(j)
  140 w(l3)=w(l3)-a(j)*w(l)
  150 continue
  160 w(l3)=w(l3)/pivot(l3)
c
c
c   compute the one-norm of vector   z   (vector   z   is
c   the vector calculated above and stored in array   w .
c
c
      znorm=0.0
      do   170   i=1,n
  170 znorm=znorm+abs(w(i))
c
c
c   find the value of the required estimate for the reciprocal
c   of the condition number of matrix   a .
c
c
      rcond=(ynorm/anorm)/znorm
c
c
c   end of the computations.
c
c
  180 if(ifail.ne.0) rcond=-1.0
      return
      end
