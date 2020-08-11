C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
