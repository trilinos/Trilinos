C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
      subroutine centroid(x, y , z, xcent, ycent, zcent,
     *  nelblk, numelb, numlnk, link, ndim)
      real x(*), y(*), z(*)
      real xcent(*), ycent(*), zcent(*)
      integer numelb(*), numlnk(*), link(*)

      IELNK = 0
      IE = 0
      DO 100 IELB = 1, nelblk
         IS = IE + 1
         IE = IE + NUMELB(IELB)
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)

         CALL cent1(x, y, z, xcent(is), ycent(is), zcent(is),
     *     NUMELB(IELB), NUMLNK(IELB), LINK(ISLNK), NDIM)
  100 CONTINUE

      RETURN
      END

      subroutine cent1(x, y, z, xcent, ycent, zcent,
     *  numelb, numlnk, link, ndim)

      real x(*), y(*), z(*)
      real xcent(*), ycent(*), zcent(*)
      integer numelb, numlnk
      integer link(numlnk,*)
      integer ndim

      do 20 ne=1, numelb
        xc = 0.0
        yc = 0.0
        zc = 0.0
        do 10 j=1, numlnk
          xc = xc + x(link(j,ne))
          yc = yc + y(link(j,ne))
          if (ndim .eq. 3) then
            zc = zc + z(link(j,ne))
          end if
 10     continue

        rnodes = numlnk
        xcent(ne) = xc / rnodes
        ycent(ne) = yc / rnodes
        if (ndim .eq. 3) then
          zcent(ne) = zc / rnodes
        end if

 20   continue
      return
      end
