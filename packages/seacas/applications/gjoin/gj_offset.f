C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      subroutine offset (off, scale, crd, length)
      real crd(length)

      do 10 i=1, length
         crd(i) = scale * crd(i) + off
 10   continue

      return
      end
