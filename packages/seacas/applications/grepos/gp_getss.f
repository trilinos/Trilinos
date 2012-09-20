C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

      subroutine getss(ndbin, numess, idess, neess, ndess, ixeess, 
     *  ixdess, lteess, ltsess, ltssnc, kfacss, a, allone, lessdf,
     *  *)
      
      
      integer idess(*), neess(*), ndess(*), ixeess(*), ixdess(*)
      integer lteess(*), ltsess(*)
      integer ltssnc(*)
      real    a(*)
      logical allone

      character*256 STRING
      
C ... Get the sideset ids
      call exgssi(ndbin, idess, ierr)
      if (ierr .gt. 0) return 1
      
      isoff = 1
      inoff = 1
      idoff = 1
      allone = .TRUE.

      call exgcssc(ndbin, ltssnc, ierr)
      if (ierr .gt. 0) then
        call exopts (EXVRBS, ierr1)
        call exerr('grepos/getss', 'Error from exgssc', ierr)
        return 1
      end if

C ... Read each sideset      
      itotdf = 0
      do 104 i=1, numess

C ... Get the sideset parameters
        call exgsp(ndbin, idess(i), nsides, ndist, ierr)
        if (ierr .gt. 0) then
           call exopts (EXVRBS, ierr1)
           call exerr('grepos/getss', 'Error from exgsp', ierr)
           return 1
        end if
        
        if (nsides .gt. 0) then
          call exgss(ndbin, idess(i), lteess(isoff), ltsess(isoff),
     *      ierr)
          if (ierr .gt. 0) then
            call exopts (EXVRBS, ierr1)
            call exerr('grepos/getss', 'Error from exgss', ierr)
            return 1
          end if
        end if
        
C ... Count number of nodes in the node list
        incnt = 0
        do 102 ii = 0, nsides-1
          incnt = incnt + ltssnc(isoff+ii)
 102    continue

        ixeess(i) = isoff
        neess(i) = nsides
        ixdess(i) = inoff
        ndess(i) = incnt
        itotdf = itotdf + incnt
        isoff  = isoff + nsides
        inoff  = inoff + incnt
 104  continue

      call mdlong('FACESS', kfacss, itotdf)
      lessdf = itotdf
      
      inoff = 0
      do i=1, numess
         call exgsp(ndbin, idess(i), nsides, ndist, ierr)
         
C     ... See if the number of distribution factors is 0 or 'incnt'
         if (ndist .eq. 0 .or. ndist .ne. ndess(i)) then
C     ... Fill in the missing dist factors with '1.0'
            do ii = 0, ndess(i)-1
               a(kfacss+inoff+ii) = 1.0
            end do
         else 
C ... Read in the existing distribution factors
            call exgssd(ndbin, idess(i), a(kfacss+inoff), ierr)
         end if
C ... There are distribution factors, but an incorrect number
         if (ndist .gt. 0 .and. ndist .ne. ndess(i)) then
            write (string, 900) idess(i)
            call sqzstr(string, lstr)
            call prterr('ERROR',STRING(:LSTR))
            write (string, 901) ndist, ndess(i)
            call sqzstr(string, lstr)
            call prterr('ERROR',STRING(:LSTR))
         end if
        
        inoff = inoff + ndess(i)
      end do
      
C... Check whether all distribution factors are constant for each sideset.
      lessnl = inoff - 1
      if (lessnl .gt. 0) then
        value = a(kfacss)
      else
        value = 1.0
      end if
      do 105 i=1, lessnl-1
         if (a(kfacss+i) .ne. value) then
            allone = .FALSE.
            go to 106
         end if
 105  continue
 106  continue
      return

 900  format('Sideset ',i5,' has distribution factor/sideset',
     $     ' node count mismatch.')
 901  format('Database has ',i7,' df, but ',i7,
     *  ' sideset nodes. All DF set to 1.0')
      end
