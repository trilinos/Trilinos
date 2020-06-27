C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: instr.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:43  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              instr
c
c DESCRIPTION:          Returns .TRUE. if a command file was specified
c                       on the BLOT command line, indicating
c                       that instructions are to be read from that file.
c                       Returns .FALSE. if no file was specified on the
c                       BLOT command line, indicating
c                       that instructions are to be read interactively.
c
c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   logical function
c
c USAGE:               instr()
c
c PARAMETERS:           none
c
c CALLS:                exname ( SUPES )
c
c GLOBAL VARIABLES REFERENCED:   none
c
c CALLING ROUTINE(S):   getins (BLOT)
c
c SYSTEM DEPENDENCIES:  none
c
c ======================================================================
c ======================================================================
c
      logical function instr()

c
c **********************************************************************
c
c        declarations
c
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
c
c **********************************************************************
c
c        data statements
c
      data symnum / 4 /
c **********************************************************************
c **********************************************************************
c

c
c           get symbol from operating system
c
      call exname ( -symnum, name, ln )

c
c           print symbol for debugging purposes
c
c     print *, 'name = ',name,' ln = ',ln
c
      if ( name(1:ln) .eq. 'YES' ) then
         instr = .TRUE.
      else
         instr = .FALSE.
      endif

      return
      end
