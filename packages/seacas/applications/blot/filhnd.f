C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================

c ROUTINE:              filhnd

c DESCRIPTION:          Opens and closes files.

c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511

c DATE:                 December 20, 1988

c TYPE OF SUBPROGRAM:   subroutine

c USAGE:               call filhnd (unit , fil1, ecodei, ecodeo,
c                                    type, fform, facces, frecl, *)

c PARAMETERS:

c        integer unit   -- (INPUT)
c                       If > 0, specifies the logical unit to be
c                          opened.
c                       If < 0, -unit specifies the logical unit to
c                          close.
c                       If = 0, all open logical units are to be
c                          closed.

c        character type -- (INPUT)
c                       'I' if input file (status = 'old')
c                       'O' if output file (status = 'new')
c                       'U' if unknown file type (status = 'unknown')
c                       'S' if scratch file (status = 'scratch')

c        character fform -- (INPUT)
c                       'F' if formatted file
c                       'U' if unformatted file

c CALLS:

c        prterr (BLOT) --     Prints an error message if one occurred
c                             during the execution of filhnd.
c        exname (SUPES)    -- Gets the filename associated with a unit
c                             number.
c        lenstr (strlib) --   Gets the length of a string (excluding
c                             trailing blanks).

c GLOBAL VARIABLES REFERENCED:

c CALLING ROUTINE(S):         getins (BLOT)

c SYSTEM DEPENDENCIES:        none

c ======================================================================
c ======================================================================

      subroutine filhnd (unit, filn, ecodei, ecodeo, type, fform,
     &   facces, frecl, *)

c        parameters

      integer unit
c          if > 0, the logical unit of the file to open.
c          if < 0, the logical unit of the file to close.
c          if = 0, close all files
      logical ecodei
c          On input, .TRUE., if program should abort
c          if an error on opening the file occurred.
c          .FALSE. otherwise.
      logical ecodeo
c          On output, .TRUE. if the file opening occurred successfully.
c          .FALSE. otherwise.
      integer frecl
c          If a file is to be opened and it is a direct access file,
c          the file's record length.
c          If a single file ( not all files ) is to be closed,
c          this parameter  = 1 to delete the file on closing.
c          Otherwise, the default status parameter is used on the close
c          call.
      character type
c          'o' if an old file,
c          'n' if a new file,
c          's' if a scratch file
c          'u' if unknown
      character fform
c          'u' if unformatted
c          'f' if formatted
      character facces
c          'd' if direct
c           s' if sequential
      character * (*) filn
c           Name of the file to open.  If ! = ' ', then filhnd calls
c           the SUPES routine EXNAME to get the filename associated
c           with the specified unit number.

c           if unit <= 0, then all other parameters are ignored.

c        declarations

      character*2048 filnam
c           filename associated with unit
      integer lname
c           length of string filnam (less trailing blanks)
      character*11 form
c           value of form keyword in open statement
      character*10 access
c           value of access keyword in open statement
      character*7 status
c           value of status keyword in open statement
      integer lform,lstat
c           length of form and status keywords  (less trailing blanks)
      integer lacces
c           length of access keyword (less trailing blanks)
      integer ios
c           integer indicating error status of the open call.
c           = 0 no error,
c           = error code otherwise.
      character*1 cparm
c           dummy argument for call to exparm
      character tform, ttype, tacces
c           Temporary variables for storing modified values of fform,
c           type, and facces

c *****************************************************************
c *****************************************************************

c        static declarations

      logical first
      save first

      integer maxunits
      parameter (maxunits=25)
c           maximum number of units that may be open at once.
      integer numopn
c            number of currently open files
      integer opnlst(maxunits)
c            array of currently open unit
c            numbers

      save numopn, opnlst

      data first / .TRUE. /

c *********************************************************************
c *********************************************************************

      if (first) then
        numopn = 0
        first = .FALSE.
      end if

c      print *, 'numopen = ',numopn
c      if ( numopn .gt. 0 )
c     &    print *, 'list is ',(opnlst(i),i=1,numopn)

      if ( unit .gt. 0 ) then

c           open file associated with unit

c                    set open keywords

         cparm = fform
         call upcase_bl ( cparm )
         tform = cparm(1:1)
         if ( tform .eq. 'U' ) then
            form = 'unformatted'
         else if (tform .eq. 'F' ) then
            form = 'formatted'
         else
            call prterr ('PROGRAM',
     &         'Bad value for fform parameter in filhnd routine.')
            return 1
         endif
         lform = lenstr ( form )

         cparm = type
         call upcase_bl ( cparm )
         ttype = cparm(1:1)
         if ( ttype .eq. 'O' ) then
            status = 'old'
         else if ( ttype .eq. 'N' ) then
            status = 'new'
         else if ( ttype .eq. 'U' ) then
            status = 'unknown'
         else if ( ttype .eq. 'S' ) then
            status = 'scratch'
         else
            call prterr ('PROGRAM',
     &         'Bad value for type parameter in filhnd.')
            return 1
         endif
         lstat = lenstr ( status )

         cparm = facces
         call upcase_bl ( cparm )
         tacces = cparm(1:1)
         if ( tacces .eq. 'D' ) then
            access = 'direct'
         else if ( tacces .eq. 'S' ) then
            access = 'sequential'
         else
            call prterr ('PROGRAM',
     &         'Bad value for access parameter in filhnd.')
            return 1
         endif
         lacces = lenstr ( access )

c                 open file

         if ( status .ne. 'scratch' ) then

c                    get file associated with unit

            filnam = filn
            call pack ( filnam, lname )
            if ( lname .eq. 0 ) then
               call exname ( unit, filnam, lname )
            endif

            if ( access .eq. 'direct' ) then
               open ( unit=unit, file=filnam(:lname),
     &            form=form(:lform),
     &            status=status(:lstat), access=access(:lacces),
     &            recl=frecl, iostat=ios )
            else
c               print *,'filename=',filnam(:lname),'='
c               print *,'form=',form(:lform),'='
c               print *,'status=',status(:lstat),'='
               open ( unit=unit, file=filnam(:lname),
     &            form=form(:lform),
     &            status=status(:lstat), iostat=ios)
            endif

            if ( ios .ne. 0 ) then
               if ( ecodei ) then
                  call prterr ('FATAL',
     &               'Error opening file in filhnd')
                  return 1
               else
                  call prterr ('ERROR',
     &               'Error opening file in filhnd')
                  ecodeo = .FALSE.
               endif
            else
               ecodeo = .TRUE.
            endif

         else

            if ( access .eq. 'direct' ) then
               open ( unit=unit, form=form(:lform),
     &            status=status(:lstat), access=access(:lacces),
     &            recl=frecl, iostat=ios )
            else
               open ( unit=unit, form=form(:lform),
     &            status=status(:lstat), iostat=ios)
            endif

            if ( ios .ne. 0 ) then
               if ( ecodei ) then
                  call prterr ('FATAL',
     &               'Error opening file in filhnd')
                  return 1
               else
                  call prterr ('ERROR',
     &               'Error opening file in filhnd')
                  ecodeo = .FALSE.
               endif
            else
               ecodeo = .TRUE.
            endif

         endif

c                 update list of open files

         if ( ecodeo ) then
            numopn = numopn + 1
            opnlst(numopn) = unit
         endif

      else if ( unit .lt. 0 ) then

c           close file

         unit = -unit
         if ( frecl .eq. 1 ) then
            close ( unit=unit, status='delete', iostat=ios )
         else
            close ( unit=unit, iostat=ios )
         endif
         if ( ios .ne. 0 ) then
            call exname ( unit, filnam, lname )
            if ( ecodei ) then
               call prterr ('PROGRAM',
     &            'Error closing file in filhnd')
               return 1
            else
               call prterr ('PROGRAM',
     &            'Error closing file in filhnd')
               ecodeo = .FALSE.
            endif
         else
            ecodeo = .TRUE.
         endif

c           update list of open files

         if ( ecodeo ) then
            i = 1
  100       continue
            if ( (i .le. numopn) .and. (unit .ne. opnlst(i)) ) then
               i = i + 1
               go to 100
            endif

            if ( i .gt. numopn ) then
               call prterr ('PROGRAM',
     &'Closed a file in filhnd that was not on the list of open files')
               return 1
            else
               numopn = numopn - 1
               do 110 j = i, numopn
                  opnlst(j) = opnlst(j+1)
  110          continue
            endif

         endif
      else

c           close all open files

         ecodeo = .TRUE.
         do 120 i = 1, numopn

            close ( unit=opnlst(i), iostat=ios )
            if ( ios .ne. 0 ) then
               call exname ( opnlst(i), filnam, lname )
               if ( ecodei ) then
                  call prterr ('PROGRAM',
     &               'Error closing file in filhnd')
                  return 1
               else
                  call prterr ('PROGRAM',
     &               'Error closing file in filhnd')
                  ecodeo = .FALSE.
               endif
            endif
  120    continue

      endif

c      print *, 'about to exit filhnd'
c      print *, 'numopen = ',numopn
c      if ( numopn .gt. 0 )
c     &    print *, 'list is ',(opnlst(i),i=1,numopn)

      return
      end
