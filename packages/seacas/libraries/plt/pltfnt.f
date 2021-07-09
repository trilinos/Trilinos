C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTFNT(FILENM)
      CHARACTER*(*) FILENM
      CHARACTER*132 TLINE

      if (FILENM .eq. 'STKFNT') then
        call plt_stick()
      else if (FILENM .eq. 'SSRFNT') then
        call plt_sanserif()
      else if (FILENM .eq. 'ROMFNT') then
        call plt_roman()
      else
        CALL PLTFLU
        TLINE = 'Error: Unrecognized font specification '//
     *    FILENM
        CALL SIORPT('PLTFNT',TLINE,2)
      end if
      RETURN

      END
