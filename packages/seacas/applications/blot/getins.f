C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================

c ROUTINE:              getins

c DESCRIPTION:          This routine is in charge of getting user
c                       input and keeping track of where
c                       it is coming from.

c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511

c DATE:                 December 20, 1988

c TYPE OF SUBPROGRAM:   subroutine

c USAGE:               call getins (id, maxfld, nfield,
c                                    kvalue, cvalue, ivalue, rvalue,
c                                    line, iostat, prompt,
c                                    lprom, *)

c PARAMETERS:

c     character*(*) id  -- (input)
c                        = 'parse' if the input line should be read
c                        and parsed by the free field reader.  If so,
c                        input information is returned in the nfield,
c                        kvalue, cvalue, ivalue, and rvalue parameters.
c                        = 'line' if the input line is to be read as
c                        one line.  If so, the input line is returned
c                        in the parameter line.
c     integer maxfld  -- (input)
c                        Maximum number of fields returned by the
c                        free field reader (frefld).
c     integer nfield  -- (output)
c                        number of data fields returned by
c                        the free field reader (frefld).
c     integer kvalue( maxfld ) -- ( output )
c                        Translation states of the data fields
c                        returned by the free field reader.  The
c                        value of each element of this array is
c                        interpreted as follows:
c                           -1 = This is a null field.
c                            0 = This is a nonnumeric field; only
c                                cvalue contains a specified value.
c                            1 = This is a real numeric field;
c                                cvalue and rvalue contain specified
c                                values.
c                            2 = This is an integer numeric field;
c                                cvalue, rvalue, and ivalue contain
c                                specified values.
c     character*(*) cvalue( maxfld ) -- ( output )
c                         Character values of the data fields returned
c                         by the free field reader (frefld).  The data
c                         will be left justified and either
c                         blank-filled or truncated.  The value in
c                         this array is set by frefld to be blank
c                         for a null field.
c     integer ivalue( maxfld ) -- ( output )
c                         Integer values of the data fields.  The
c                         value in this array is set by frefld to
c                         be zero for a null or non-integer field.
c     real rvalue( maxfld ) -- ( output )
c                         Floating-point values of the data fields.
c                         The value in this array is set by frefld
c                         to be zero for a null or non-real fields.
c     character*(*) line -- ( output )
c                         If id = 'line', then the input line
c                         is returned in this parameter.
c     integer iostat -- (output)
c                         iostat value returned by the frefld reader
c                         (if id = 'parse') or by the read
c                         statement (if id = 'line').
c     character*(*) prompt -- (input)
c                         Prompt to be displayed to user.
c     integer lprom -- (input)
c                         Length of prompt string.
c     * -- Alternate return in case of fatal error.

c CALLS:                prterr (etclib), filhnd (BLOT),
c                       lenstr (strlib),
c                       frefld (SUPES), exname (SUPES)

c GLOBAL VARIABLES REFERENCED:

c     none

c CALLING ROUTINE(S):   comand (BLOT)

c SYSTEM DEPENDENCIES:  none

c ======================================================================
c ======================================================================

      subroutine getins (id, maxfld, nfield, kvalue, cvalue,
     &   ivalue, rvalue, line, iostat, prompt, lprom, *)

c ***********************************************************************

c        parameters

      character*(*) id
      integer maxfld
      character*(*) cvalue( maxfld )
      integer nfield, kvalue( maxfld ), ivalue( maxfld )
      real rvalue( maxfld )
      character*(*) line
      integer iostat
      character*(*) prompt
      integer lprom

c ***********************************************************************

c        local declarations

      integer maxstk
      parameter ( maxstk = 10 )
c         the maximum number of files allowed in the stack of files
c         containing instructions that are still to be read.
      character*2048 file, filelc
c         Name of file to be opened for reading instructions
      character*2048 name( maxstk )
c         the stack of file names containing instructions that are
c         still to be read
      integer recred( maxstk )
c         recred(i) = the number of records that have already been
c         processed in file i of the stack.  That is, when file i
c         is reopened as the source of instructions, the first
c         recred(i) records (instructions) have already been read
c         and so should be skipped.
      integer stkpnt
c         the number of files currently in the stack.  That is,
c         when the last instruction is read from a file,
c         the file from the top of the stack, file( numstk ),
c         is reopened.  This is the file that called the last file.
      save stkpnt, name, recred

      logical first
c         .TRUE. until after the first call to getins.
c         .FALSE. thereafter.
      save first
      integer ln
c          length of a file name or a command string
      integer ios
c          I/O status returned by and I/O statement or the free
c          field reader (frefld).
      logical gotins
c          .TRUE. if an instruction was successfully read from
c          the input stream.
c         .FALSE. otherwise. ( e.g., an EOF mark was read from
c          the input file.
      logical ecode
c          code returned by subroutine filhnd.
c          .TRUE. if file handling proceeded properly.
c          .FALSE. otherwise
      integer temp
c          temporary variable for saving integer variables.
      logical try
c          .TRUE. if an attempt to open an instruction file is to be made.
c          .FALSE. otherwise.
      logical quit
c          .TRUE. if control should go back to the top of the routine
c          to try to read another command line.
c          .FALSE. if a command line has successfully been read and
c          the program should exit.
      character*2048 cmdfile
c          Name of the command for directing command to a new file.
      save cmdfile

      logical parse
c          .TRUE. if the input is to be read and parsed by the free
c          field reader.
c          .FALSE. if the input is to be read as a line.

      logical instr

      integer nin
c          logical unit where instructions are currently being read from.
      save nin

      character*2048 cval2(80)
      logical batch

c ***********************************************************************

c     data statements

      data recred / maxstk*0 /
      data first / .TRUE. /
      data cmdfile(1:7) / 'CMDFILE' /

c ***********************************************************************
c ***********************************************************************

      iostat = 0
      if ( id .eq. 'parse' ) then
         parse = .TRUE.
      else if ( id .eq. 'line' ) then
         parse = .FALSE.
      else
         call prterr ('PROGRAM',
     &      'illegal value for id parameter passed to getins')
         return 1
      endif

      if ( first ) then
        stkpnt = 0
c         first call to getins

c           find out if an instruction file was specified in the P4
c           field of the command line

         if ( instr() ) then
c   an instruction file was specified
            call exname ( 7, name( 1 ), ln )
c   get name of file

c                 open file

            if ( batch() ) then
               call filhnd ( 7, name(1)(:ln), .TRUE., ecode,
     &            'o', 'f', 's', 0, *150)
               stkpnt = stkpnt + 1
c   increment stack pointer
               nin = 7
c    set logical unit where instructions are to be read from
            else
               call filhnd ( 7, name(1)(:ln), .FALSE., ecode,
     &            'o', 'f', 's', 0, *150)
               if ( ecode ) then
c  file opened properly
                  stkpnt = stkpnt + 1
c  increment stack pointer
                  nin = 7
c   set logical unit where instructions are to be read from
               else
c   instruction file wasn't there. prompt user for instructions.
                  nin = 0
               endif
            endif

         else
c   no instruction file was specified
            if ( batch() ) then
               call prterr ('FATAL',
     &        'No instruction file was specified on the command line.')
               return 1
            else
               nin = 0
            endif
         endif
         first = .FALSE.

      endif

c              get an instruction

  100 continue

      gotins = .FALSE.

  110 continue

      if ( nin .eq. 7 )
     &   recred( stkpnt ) = recred( stkpnt ) + 1

c              read instruction

      if ( parse ) then
         call frefld ( nin, 0, prompt(:lprom), maxfld, ios,
     &      nfield, kvalue, cval2, ivalue, rvalue )
      else
        if (nin .eq. 0) then
          read (*, fmt=10010, iostat=ios ) line
        else
          read ( unit=nin, fmt=10010, iostat=ios ) line
        end if
10010     format ( a )
      endif

c           check for an error in reading instruction

      if ( ios .gt. 0 ) then
c   error in reading instruction
         call prterr ('FATAL',
     &      'error reading an instruction in getins')
         return 1
      else if ( ios .lt. 0 ) then

c              end of file read from input stream.

         if ( nin .ne. 0 ) then
c            EOF mark read from a file.
c            Close the file and open the
c            previous one on the stack, if
c            one exists.  If not, switch to
c            interactive input.
            call filhnd ( -nin, ' ', .TRUE., ecode,
     &         ' ', ' ', ' ', 0, *150)
            stkpnt = stkpnt - 1
            if ( stkpnt .gt. 0 ) then

c                    open previous file

               call filhnd ( nin, name( stkpnt ), .TRUE., ecode,
     &            'o', 'f', 's', 0, *150)

c                    skip records that have already been read

               do 120 i = 1, recred( stkpnt )
                  read ( nin, * )
  120          continue

            else
c  no more files on the stack.  Switch
c  to interactive input if an interactive
c  job.  If a batch job, print error.
               if ( .not. batch() ) then
                  nin = 0
               else
                  call prterr ('PROGRAM',
     &          ' End of file mark read in lowest level command file.')
                  return 1
               endif

            endif

         endif

      else
         gotins = .TRUE.
      endif

      if ( .not. gotins ) go to 110

c              if id = 'line', just return the line of input

      if ( .not. parse ) then

         try = .FALSE.
         quit = .TRUE.

      else

c              check if command specifies that instructions
c              are to be read from a different source

         ln = lenstr ( cval2(1) )
         if ( cval2(1)(1:ln) .eq. cmdfile(1:ln) ) then

            if ( kvalue(2) .eq. 0 ) then
               file = cval2(2)
C ... Convert filename to all lowercase -- FREFLD converts to all uppercase
               call lowstr(filelc, file)
               try = .TRUE.
            else
               if ( nin .eq. 0 ) then
                  write ( *, 10020 )
10020              format ( // 5x, 'Please enter an instruction file',
     &               ' to read from or ''a'' to ' / 5x,
     &               'continue reading from the current unit. > ', $ )
                  read ( *, 10030 ) file
10030              format ( a )
               else
                  recred( stkpnt ) = recred( stkpnt ) + 1
                  read ( nin, 10030 ) file
               endif
               call pack ( file, ln )
               if ( file(1:1) .eq. 'A'  .or.
     &            file(1:1) .eq. 'a' ) then
                  try = .FALSE.
                  quit = .FALSE.
               else
                  try = .TRUE.
               endif
            endif
         else
            try = .FALSE.
            quit = .TRUE.
            do 130 i2 = 1, maxfld
               cvalue(i2) = cval2(i2)
  130       continue
         endif

      endif

      if ( try ) then

         quit = .FALSE.
c          will want to return to the
c          top of the routine to read
c          an instruction.

         temp = nin
c          remember where instructions are coming from.

         if ( nin .eq. 7 ) then
c           instructions are currently coming
c           from a file.  This current file
c           must be closed.  The information
c           on this file ( its name and the
c           pointer to its current instruction
c           has already been stored on the stack.
            call filhnd ( -nin, ' ', .TRUE., ecode, ' ', ' ', ' ',
     &         0, *150)
         else
c           instructions are currently coming
c           from the terminal.  Change unit
c           of instructions to 7
            nin = 7
         endif

         stkpnt = stkpnt + 1
c            increment stack pointer

         if ( stkpnt .gt. maxstk ) then
            call prterr ('FATAL',
     &'Nesting of instruction files is greater than maximum allowed')
            return 1
         endif

c              open instruction file

C ... We have a problem on systems with case-sensitive file names.
C     FREFLD converts all strings to uppercase.  Therefore, the filename
C     specified in 'file' will be all uppercase.  Since the file is
C     typically lowercase, we lowercase 'file' in variable 'filelc'.
C     First, try to open 'file' (uppercase), if this fails, try to
C     open 'filelc' (lowercase).  NOTE: There is no way (yet) to deal
C     with mixed-case filenames.  GDS 7/1/91.

         if ( batch() ) then
            call filhnd ( nin, file, .FALSE., ecode,
     &         'o', 'f', 's', 0, *150)
            if (.not. ecode) then
               call filhnd ( nin, filelc, .TRUE., ecode,
     &              'o', 'f', 's', 0, *150)
            end if
         else
            call filhnd ( nin, file, .FALSE., ecode,
     &         'o', 'f', 's', 0, *150)

            if (.not. ecode) then
               call filhnd ( nin, filelc, .FALSE., ecode,
     &              'o', 'f', 's', 0, *150)
            end if
         endif
         if ( ecode ) then
c   file was successfully opened
            name( stkpnt ) = file
            recred( stkpnt ) = 0
         else
c   file was not successfully opened
c   reopen last instruction file ( if
c   necessary ) and reset parameters.
            call prterr ('ERROR',
     &         ' Command file could not be opened.  Command ignored. ')

            nin = temp
            stkpnt = stkpnt - 1
            if ( nin .eq. 7 ) then

               call filhnd ( nin, name( stkpnt ), .TRUE., ecode,
     &            'o', 'f', 's', 0, *150)

c                       skip records that have already been read

               do 140 i = 1, recred( stkpnt )
                  call frefld ( nin, 0, ' ', maxfld, ios, nfield,
     &               kvalue, cval2, ivalue, rvalue )
                  if ( ios .ne. 0 ) then
                     call prterr ('PROGRAM',
     & 'Error skipping previously read records in getins')
                     return 1
                  endif
  140          continue
            endif
         endif

      endif

c              return to the top of the routine to read
c              another instruction, if necessary

      if ( .not. quit ) go to 100
      return

c              alternate return

  150 continue
      return 1
      end
