C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C************************************************************************
      subroutine version (qainfo)
C************************************************************************

      include 'exodusII.inc'
      character*(mxstln) qainfo(6)

C      --QAINFO - the current program QA information:
C      --   (1) = program name
C      --   (2) = revision date
C      --   (3) = version as "QA xx.xx" or "X  xx.xx" or "   xx.xx"
C      --   (4) = program name with version appended
C      --   (5) = date of current run
C      --   (6) = time of current run

      qainfo(1) = 'mapvar                          '
      qainfo(2) = '2019/05/15                      '
      qainfo(3) = ' 2.01                           '
      qainfo(4) = '                                '
      qainfo(5) = '                                '
      qainfo(6) = '                                '
      call exdate(qainfo(5))
      call extime(qainfo(6))

      return
      end
