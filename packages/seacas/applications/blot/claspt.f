C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: claspt.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1996/06/21 16:07:01  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:55:47  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:10  gdsjaar
c Added RCS Id and Log to all files
c
      subroutine claspt( xpt, ypt, zpt, cutpt, cutnrm, status)

      real xpt, ypt, zpt
      real cutpt(3), cutnrm(3)
      integer status
      parameter(ISIN=1, ISON=-2, ISOUT=-1)
      real vec(3)
      real tol
      parameter(REFTOL=1e-4)
c
c check dot product of normal vector and (pt-cutpt) vector to find
c if point is in front or behind plane
c
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
c
      return
      end
