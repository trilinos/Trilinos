      subroutine   corsum   ( tmmtch, timhrz, timvrt, nhrows, nhcols,
     $                        nsrows, nscols, nvrows, nvcols, output )

c     ==================================================================
c     ==================================================================
c     ====  corsum -- print summary of coarse decomposition         ====
c     ==================================================================
c     ==================================================================

c     created by john lewis, boeing computer services, sept. 17, 1990

      integer            nhcols, nhrows, nscols, nsrows, nvcols, nvrows,
     $                   output

      real               tmmtch, timhrz, timvrt

c     ==================================================================

      write (output, 60000) tmmtch, timhrz, timvrt, nhrows, nhcols,
     $                      nsrows, nscols, nvrows, nvcols

60000 format (/ '0coarse decomposition',
     $        / '0   time required for:',
     $        / '       maximum matching:', 1pe10.1,
     $        / '       horizontal block:', 1pe10.1,
     $        / '       vertical block  :', 1pe10.1,
     $        / '0   block dimensions --  rows   columns',
     $        / '         horizontal', 2i9, 
     $        / '         square    ', 2i9, 
     $        / '         vertical  ', 2i9 )

      return
      end

