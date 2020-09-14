C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
      subroutine elementize(var_nod, var_el,
     *  nelblk, numelb, numlnk, link)
      real var_nod(*), var_el(*)
      integer numelb(*), numlnk(*), link(*)

      IELNK = 0
      IE = 0
      DO 100 IELB = 1, nelblk
         IS = IE + 1
         IE = IE + NUMELB(IELB)
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)

         CALL elemtz1(var_nod, var_el(is),
     *     NUMELB(IELB), NUMLNK(IELB), LINK(ISLNK))
  100 CONTINUE

      RETURN
      END

      subroutine elemtz1(var_nod,  var_el, numelb, numlnk, link)

      real var_nod(*)
      real var_el(*)
      integer numelb, numlnk
      integer link(numlnk,*)

      do 20 ne=1, numelb
        var = 0.0
        do 10 j=1, numlnk
          var = var + var_nod(link(j,ne))
 10     continue

        rnodes = numlnk
        var_el(ne) = var / rnodes
 20   continue
      return
      end
