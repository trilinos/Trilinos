C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine mapvar(nold, nnew, nvar, map, vars, scr)
C=======================================================================
C   --*** MAPVAR *** (GREPOS) Map values from old array to new.
C   --
C   --MAPVAR reorders the data in VARS based on the map in MAP
C   --
C   --Parameters:
C   --   NOLD - IN - number of old values,
C   --   NNEW - IN - number of new values
C   --   NVAR - IN - the number of variables
C   --   NVARDM - IN - the row dimension of VARNP
C   --   MAP  - IN - map from new value to old MAP(NEW) = OLD
C                    size is 'nnew'
C   --   VARS - IN/OUT - the values. On input in old order,
C                        on output in new order
C   --   SCR  - TMP - temporary storage area

      integer map(*)
      real vars(*)
c     real vars(nold, nvar)
      real scr(*)

C ... VARS should really be a doubly-dimensioned array (NOLD, NVAR),
C     The dimensions need to be in this order so we can read them
C     in using exgev and exgnv.  But, this order doesn't work very
C     well when trying to reduce the size of NOLD

C ... TODO: Need to use the truth table to make sure variables
C           exist for each element.
      do 30 ivar = 1, nvar
        do 10 i = 1, nnew
          scr(i) = vars(map(i) + nold * (ivar-1) )
 10     continue

        do 20 i = 1, nnew
          vars(i + nnew * (ivar-1)) = scr(i)
 20     continue
 30   continue

      return
      end
