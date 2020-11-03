C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      subroutine claspt( xpt, ypt, zpt, cutpt, cutnrm, status)

      real xpt, ypt, zpt
      real cutpt(3), cutnrm(3)
      integer status
      parameter(ISIN=1, ISON=-2, ISOUT=-1)
      real vec(3)
      real tol
      parameter(REFTOL=1e-4)
c check dot product of normal vector and (pt-cutpt) vector to find
c if point is in front or behind plane
      vec(1) = xpt - cutpt(1)
      vec(2) = ypt - cutpt(2)
      vec(3) = zpt - cutpt(3)
c      tol = amax1( vec(1), vec(2), vec(3) ) * REFTOL
      tol = REFTOL
      dot = vec(1)*cutnrm(1) + vec(2)*cutnrm(2) + vec(3)*cutnrm(3)

      if( abs(dot) .lt. tol) then
          status = ISON
      else if(dot .gt. 0.0) then
          status = ISOUT
      else
          status = ISIN
      end if
      return
      end
