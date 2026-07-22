C Copyright(C) 1999-2020, 2024, 2025 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE VERSION(QAINFO)

      include 'params.blk'

      CHARACTER*(MXQARC) QAINFO(6)

      QAINFO(1) = 'blot                            '
      QAINFO(2) = '2025/01/03                      '
      QAINFO(3) = ' 3.14159                        '
      QAINFO(4) = '                                '

      RETURN
      END
