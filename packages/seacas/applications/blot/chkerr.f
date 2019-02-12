C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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

C $Log: chkerr.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1997/11/11 14:55:53  gdsjaar
C Added 'external blkdat' to main program to ensure that the block data
C gets linked into the executable. Wasn't happening on dec alpha
C systems.
C
C Removed unreachable lines in several routines
C
C Fixed variable name spelling in contor.f
C
C Unsplit strings that were split across lines
C
C Removed old error variables left over from exodusIIv1
C
C Upped version number
C
C Revision 1.1  1994/04/07 19:55:29  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/07/28  15:58:32  gdsjaar
c Fixed junky looking output
c
c Revision 1.2  1993/07/27  20:27:18  gdsjaar
c Added error checking, set to non-verbose exodus II file opening,
c cleaned up dbiv0 routine, removed some unused variables.
c
c Revision 1.1  1993/07/27  19:35:17  gdsjaar
c Added file for error checking in exodus II calls.
c
c Revision 1.3  1993/06/24  15:27:20  gdsjaar
c Added more error flags and messages.
c
c Revision 1.2  1993/03/03  17:44:44  gdsjaar
c Fixed problem with assumed length strings.
c
c Revision 1.1  1992/06/08  22:23:07  gdsjaar
c New routine to check error status of exo2 library calls
c
C=======================================================================
      subroutine chkerr (routine, caller, ierr)
C=======================================================================
      character*(*) routine, caller
      include 'exodusII.inc'
      character*80 string, path
      character*40 errors(19)
      integer ierr
      integer idebug
      data errors /
     *  'Not a netcdf id',
     *  'Too many netcdfs open',
     *  'netcdf file exists && NCNOCLOB',
     *  'Invalid Argument',
     *  'Write to read only',
     *  'Operation not allowed in data mode',
     *  'Operation not allowed in define mode',
     *  'Coordinates out of Domain',
     *  'MAXNCDIMS exceeded',
     *  'String match to name in use',
     *  'Attribute not found',
     *  'MAXNCATTRS exceeded',
     *  'Not a netcdf data type',
     *  'Invalid dimension id',
     *  'NCUNLIMITED in the wrong index',
     *  'MAXNCVARS exceeded',
     *  'Variable not found',
     *  'Action prohibited on NCGLOBAL varid',
     *  'Not a netcdf file'/

      idebug = 1
      write (path, '(A,A2,A)') CALLER, '->', ROUTINE
      lp = lenstr(path)
      if (idebug .gt. 0 .and. ierr .gt. 0)
     *  call prterr ('CMDSPEC', PATH(:LP))

      if (ierr .eq. 0) then
         continue
      else if (ierr .eq. EXBFID) then
         write (string, '(A,A)')
     $        'The EXODUS file ID is incorrect in ', path(:lp)
         call prterr ('FATAL', string)
         stop 'EXBFID'
      else if (ierr .eq. EXWARN) then
         write (string, '(A,A)')
     $        'A non-fatal error occurred in ', path(:lp)
         call prterr ('WARNING', string)
      else
         if (ierr .le. 19) then
           l = lenstr(errors(ierr))
           write (string, '(A,A,A)') errors(ierr)(:l), ' ', path(:lp)
         else
           write (string, '(A,A)')
     *       'Unknown Error ', path(:lp)
         end if
         call prterr ('WARNING', string)
      end if

      return
      end
