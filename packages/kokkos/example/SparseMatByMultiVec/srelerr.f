      double precision function srelerr ( n, xtrue, xcomp )
c
c Computes the relative error || xtrue - xcomp ||
c                             -------------------
c                                  || xtrue ||
c
c  where the norm is the infinity norm.
c  n     - Dimension
c  xtrue - Exact solution
c  xcomp - Approximate solution
c
      implicit double precision (a-h,o-z)
      integer n, i
      real*8 xtrue(n), xcomp(n), diff, tmp, diffmx, tmpmx
c
      diffmx = 0.0D0
      tmpmx = 0.0D0
      do 10 i=1,n
        diff = abs(xtrue(i)-xcomp(i))
        tmp = abs(xtrue(i))
        if (diff.gt.diffmx)  diffmx = diff
        if (tmp.gt.tmpmx)  tmpmx = tmp
 10   continue
      srelerr = diffmx/tmpmx
      return
      end
