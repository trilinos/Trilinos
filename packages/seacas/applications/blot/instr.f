C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================

c ROUTINE:              instr

c DESCRIPTION:          Returns .TRUE. if a command file was specified
c                       on the BLOT command line, indicating
c                       that instructions are to be read from that file.
c                       Returns .FALSE. if no file was specified on the
c                       BLOT command line, indicating
c                       that instructions are to be read interactively.

c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511

c DATE:                 December 20, 1988

c TYPE OF SUBPROGRAM:   logical function

c USAGE:               instr()

c PARAMETERS:           none

c CALLS:                exname ( SUPES )

c GLOBAL VARIABLES REFERENCED:   none

c CALLING ROUTINE(S):   getins (BLOT)

c SYSTEM DEPENDENCIES:  none

c ======================================================================
c ======================================================================

      logical function instr()

c **********************************************************************

c        declarations

      integer symnum
c          the number of the system command symbol defining if a
c          command file was specified.
      character*12 name
c          string that holds the value of the
c          symbol associated with the symnum
c          = "YES" if a file was specified
c            in the P4 field of the TRINITY
c            command line, indicating that
c            instructions are to be read
c            from that file
c          = "NO" otherwise.
      integer ln
c          length of the string stored in name

c **********************************************************************

c        data statements

      data symnum / 4 /
c **********************************************************************
c **********************************************************************

c           get symbol from operating system

      call exname ( -symnum, name, ln )

c           print symbol for debugging purposes

c     print *, 'name = ',name,' ln = ',ln

      if ( name(1:ln) .eq. 'YES' ) then
         instr = .TRUE.
      else
         instr = .FALSE.
      endif

      return
      end
