C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE VERSION(QAINFO)

      include 'params.blk'

      CHARACTER*(MXQARC) QAINFO(6)

      QAINFO(1) = 'blot                            '
      QAINFO(2) = '2019/03/18                      '
      QAINFO(3) = ' 3.14                           '
      QAINFO(4) = '                                '

      RETURN
      END
