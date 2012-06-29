C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Log: filhnd.f,v $
C Revision 1.4  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  2002/11/27 16:19:09  gdsjaar
C Fix filhnd calls to not pass partially uninitialized character strings to upcase.
C
C Revision 1.2  1997/11/11 14:55:55  gdsjaar
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
C Revision 1.1  1994/04/07 20:00:47  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:32  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              filhnd
c
c DESCRIPTION:          Opens and closes files.
c
c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   subroutine
c
c USEAGE:               call filhnd (unit , fil1, ecodei, ecodeo,
c                                    type, fform, facces, frecl, *)
c
c PARAMETERS:
c
c        integer unit   -- (INPUT)
c                       If > 0, specifies the logical unit to be
c                          opened.
c                       If < 0, -unit specifies the logical unit to
c                          close.
c                       If = 0, all open logical units are to be
c                          closed.
c
c        character type -- (INPUT)
c                       'I' if input file (status = 'old')
c                       'O' if output file (status = 'new')
c                       'U' if unknown file type (status = 'unknown')
c                       'S' if scratch file (status = 'scratch')
c
c        character fform -- (INPUT)
c                       'F' if formatted file
c                       'U' if unformatted file
c
c CALLS:
c
c        prterr (BLOT) --     Prints an error message if one occurred
c                             during the execution of filhnd.
c        exname (SUPES)    -- Gets the filename associated with a unit
c                             number.
c        lenstr (strlib) --   Gets the length of a string (excluding
c                             trailing blanks).
c
c GLOBAL VARIABLES REFERENCED:
c
c CALLING ROUTINE(S):         getins (BLOT)
c
c SYSTEM DEPENDENCIES:        none
c
c ======================================================================
c ======================================================================
c
      subroutine filhnd (unit, filn, ecodei, ecodeo, type, fform,
     &   facces, frecl, *)
c
c
c        parameters
c
      integer unit               
c          if > 0, the logical unit of the file to open.
c          if < 0, the logical unit of the file to close.
c          if = 0, close all files
      logical ecodei             
c          On input, .TRUE., if program should abort
c          if an error on opening the file occured.
c          .FALSE. otherwise.
      logical ecodeo             
c          On output, .TRUE. if the file opening occured successfully.
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
c
c           if unit <= 0, then all other paramaters are ignored.
c
c
c        declarations
c
      character*132 filnam       
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
c
c *****************************************************************
c *****************************************************************

c
c        static declarations
c
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
c
c           open file associated with unit
c

c                    set open keywords
c
         cparm = fform
         call upcase ( cparm )
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
         call upcase ( cparm )
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
c
c
         cparm = facces
         call upcase ( cparm )
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
c
c
c                 open file
c
         if ( status .ne. 'scratch' ) then
c
c                    get file associated with unit
c
            filnam = filn
            call pack ( filnam, lname )
            if ( lname .eq. 0 ) then
               call exname ( unit, filnam, lname )
            endif
c
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

c
c
c                 update list of open files
c
         if ( ecodeo ) then
            numopn = numopn + 1
            opnlst(numopn) = unit
         endif
c
c

      else if ( unit .lt. 0 ) then
c
c
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
c
c           update list of open files
c
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
c
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
