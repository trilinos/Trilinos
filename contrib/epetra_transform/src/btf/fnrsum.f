      subroutine   fnrsum   ( fstcmp, lstcmp, ccmstr, rcmstr, output )

c     ==================================================================
c     ==================================================================
c     ====  fnrsum -- print fine summary of fine structure          ====
c     ==================================================================
c     ==================================================================

c     created by john lewis, boeing computer services, sept. 17, 1990

      integer       fstcmp, lstcmp, ccmstr (*), rcmstr (*), output

      integer       i, cmpcol, cmpnum, cmprow

c     ==================================================================

      cmpnum = 1

      do 100 i = fstcmp, lstcmp

         cmpcol = ccmstr(i+1) - ccmstr(i)
         cmprow = rcmstr(i+1) - rcmstr(i)

         write (output, '(15x, i12, 2i10)') cmpnum, cmprow, cmpcol
         cmpnum = cmpnum+1

 100  continue

      return
      end

