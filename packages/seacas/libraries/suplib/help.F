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

C=======================================================================
      LOGICAL FUNCTION HELP (PROG, TYPE, OPTION)
C=======================================================================
C$Id: help.F,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: help.F,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1998/04/16 04:25:07  gdsjaar
CRenamed c_xhelp to cxhelp to avoid problems with f2c and g77 (add double underscore to names containing an underscore)
C
CRevision 1.1  1997/03/20 22:46:55  caforsy
Ccreated help.F from help.f
C
CRevision 1.3  1991/04/16 16:39:42  gdsjaar
CReturned to .f from help.F
C
c Revision 1.2  1991/04/08  22:06:08  gdsjaar
c Fixed bug in lenstr call
c
c Revision 1.1  1991/04/08  19:12:50  gdsjaar
c Initial revision
c
c Revision 1.2  1990/11/30  09:50:57  gdsjaar
c Modified to work on Unicos
c
c Revision 1.1.1.1  90/08/14  16:14:50  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:14:49  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:28  gdsjaar
c Initial revision
c 

C   --*** HELP *** (ETCLIB) Display HELP for a program
C   --   Written by Amy Gilkey - revised 08/27/87
C   --   D. P. Flanagan, 12/18/83
C   --
C   --HELP displays help for a program in a system dependent fashion.
C   -- !#VAX
C   --It assumes that a help command has been established in ACC$HELP !#VAX
C   --with an entry {PROG} {TYPE}. !#VAX
C   --
C   --Parameters:
C   --   HELP - IN - returns true iff the system has some form of HELP
C   --   PROG - IN - the program name
C   --   TYPE - IN - the HELP type (e.g., "COMMANDS")
C   --   OPTION - IN - the TYPE option (e.g., specific command)

C   --System Dependencies (VMS): !#VAX
C   --   LBR$OUTPUT_HELP, LIB$PUT_OUTPUT, LIB$GET_INPUT are used to !#VAX
C   --      obtain on-line help. !#VAX
C#CTSSC   --System Dependencies (CTSS):
C#CTSSC   --   No help available.

C   --Routines Called:
C   --   LBR$OUTPUT_HELP - Obtain on-line help !#VAX
C   --   LENSTR - (STRLIB) Find string length !#VAX
C   --   SQZSTR - (STRLIB) Remove extra blanks from string !#VAX

      CHARACTER*(*) PROG
      CHARACTER*(*) TYPE
      CHARACTER*(*) OPTION

      CHARACTER*80 STRING 

      HELP = .TRUE.

      IF (HELP) THEN
#ifdef USEHELP
         WRITE (*, 10)
         string = prog(:lenstr(prog)) // ' ' //
     $            type(:lenstr(type)) // ' ' //
     $          option(:lenstr(option))

         CALL CXHELP(string, lenstr(string))
#else
         WRITE (*, 10)
     &      'Help for ', PROG(:LENSTR(PROG)), ' is not available'
#endif
      ELSE
         WRITE (*, 10)
     &      'Help for ', PROG(:LENSTR(PROG)), ' is not available'
      END IF

      RETURN

   10 FORMAT (/, 1X, 5A)
      END

