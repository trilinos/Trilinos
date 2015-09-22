C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
      SUBROUTINE MDINIT (USERV)
C
C     THIS SUBPROGRAM PROVIDES DATA MANAGEMENT WITH DYNAMIC MEMORY
C     ALLOCATION CAPABILITIES.  THE ENTIRE PACKAGE IS WRITTEN IN
C     STANDARD FORTRAN-77 WITH THE EXCEPTION OF TWO SYSTEM DEPENDENT
C     SUBPROGRAMS:  IXLNUM AND EXMEMY.
C
C     AUTHOR:  W. C. MILLS-CURRAN, SANDIA NATIONAL LABORATORIES,
C     DIV. 1523.
C     DATE:  FEBRUARY, 1986.
C
C     DOCUMENTATION:  SAND86-0911
C
C     REVISED 28-OCT-1988
C     Modified to allow memory for character variables.  The error
C     handling was also improved to allow better user control.
C
C ****************************************************************************
C
C                   ISSUED BY SANDIA NATIONAL LABORATORIES,
C                          A PRIME CONTRACTOR TO THE
C                      UNITED STATES DEPARTMENT OF ENERGY
C
C ****************************************************************************
C
C THIS CODE WAS PREPARED IN THE COURSE OF WORK SPONSORED BY THE UNITED STATES
C GOVERNMENT.  NEITHER THE UNITED STATES, NOR THE UNITED STATES DEPARTMENT OF
C ENERGY, NOR THE UNITED STATES NUCLEAR REGULATORY COMMISSION, NOR ANY OF
C THEIR EMPLOYEES, NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR
C EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR USEFULNESS OF
C ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS DISCLOSED, OR REPRESENTS THAT
C ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C ****************************************************************************
C
C     ENTRY POINTS:
C
C     MDINIT (USERV)
C     MCINIT (UCHAR)
C     Initializes the data manager and calculates the offset of the
C     users base vector.
C
C     MDCOMP
C     MCCOMP
C     Compresses all storage to consolidate the voids within a
C     memory block.
C
C     MDDEBG (UNIT)
C     MCDEBG (UNIT)
C     Turns on immediate error output at the time an error is triggered.
C
C     MDDEL (NAME1)
C     MCDEL (NAME1)
C     "Deletes" (frees the storage of) an array.
C
C     MDEFIX (ENUM, COUNT)
C     MCEFIX (ENUM, COUNT)
C     Resets the error counter for a particular error flag.
C
C     MDEROR (UNIT)
C     MCEROR (UNIT)
C     Prints an error summary.
C
C     MDERPT (ENUM, COUNT)
C     MCERPT (ENUM, COUNT)
C     Returns the error count for a particular error flag.
C
C     MCEXEC ()
C     MDEXEC ()
C     Turns off deferred allocation and resolves all deferred requests.
C
C     MDFILL (TFILL)
C     MCFILL (NAME1)
C     Turns on data initialization to a specific value.
C
C     MDFIND (NAME1, NEWLOC, NEWLEN)
C     MCFIND (NAME1, NEWLOC, NEWLEN)
C     Reports the location and length of an array.
C
C     MDFOFF ()
C     MCFOFF ()
C     Turns off data initialization.
C
C     MDGET (MNGET)
C     MCGET (MNGET)
C     Requests a contiguous block of storage.
C
C     MDGIVE
C     MCGIVE
C     Returns storage to the system under certain conditions.
C
C     MDLAST (ENUM)
C     MCLAST (ENUM)
C     Returns the flag number of the last error.
C
C     MDLIST (UNIT)
C     MCLIST (UNIT)
C     Prints a storage summary.
C
C     MDLONG (NAME1, NEWLOC, NEWLEN)
C     MCLONG (NAME1, NEWLOC, NEWLEN)
C     Changes the length of an array.
C
C     MDNAME (NAME1, NAME2)
C     MCNAME (NAME1, TNAME2)
C     Changes the name of an array.
C
C     MDMEMS (NSUA, NSUD, NSUV, NSULV)
C     Reports numeric storage information.
C
C     MCMEMS (NSUA, NSUD, NSUV, NSULV)
C     Reports character storage information.
C
C     MDPRNT (NAME1, UNIT, NAME2)
C     MCPRNT (NAME1, UNIT, NEWLEN)
C     Prints an array.
C
C     MDRSRV (NAME1, NEWLOC, NEWLEN)
C     MCRSRV (NAME1, NEWLOC, NEWLEN)
C     Reserves space for an array and reports the storage location.
C
C     MDSTAT (MNERRS, MNUSED)
C     MCSTAT (MNERRS, MNUSED)
C     Reports on accumulated errors and space used.
C
C     MCWAIT ()
C     MDWAIT ()
C     Turns on deferred allocation mode.
C
C     Associated subroutines:
C
C     IXLCHR            MXLIST             MYLONG
C     MDINIT            MXLONG             MYMEMY
C     MEMTEST           MXLOOK             MYNSRT
C     MXCOMP            MXNSRT             MYPRNT
C     MXDEL             MXPRNT             MYRSRV
C     MXEROR            MXRSRV             SHFTC
C     MXEXEC            MYCOMP             SHFTI
C     MXFIND            MYDEL              SRCHC
C     MXGET             MYFIND             SRCHI
C     MXGIVE            MYGET              VTABLE
C
C
      IMPLICIT INTEGER (A-Z)
C
C     Dictionary and tables lengths.
C
C     LDICT     Maximum number of dictionary entrys/allocated arrays.
C
C     LVOID     Maximum number of voids.
C
C     These parameters may be changed independently.
C
      PARAMETER (LPDICT=256, LPVOID=2048)
      CHARACTER*8 DICT(LPDICT, 2)
      DIMENSION DPOINT(LPDICT, 2, 3)
C                              *  *> 1: LOCATION
C                              *     2: NUMERIC STORAGE UNITS
C                              *     3: CHARACTER STORAGE UNITS
C                              *        (-1 FOR NUMERIC)
C                              *> 1: NUMERIC
C                                 2: CHARACTER FOR NON MIXED MODE
      DIMENSION VOID(LPVOID, 2, 2)
C                            *  *> 1: LOCATION
C                            *     2: SIZE (IN NUMERIC UNITS FOR MIXED MODE)
C                            *> 1: NUMERIC
C                               2: CHARACTER FOR NON MIXED MODE
C
C     The strategy for the above arrays is to fix the second dimension
C     at 1 for mixed mode, but let the first dimension "overrun" so
C     that the storage is usable.  In nonmixed mode, the first dimension
C     will not overrun and the character stuff will use a value of 2 for
C     for the second dimension.  The mixed vs. nonmixed mode is
C     indicated by the variable CHRCOL. (1 = mixed, 2 = nonmixed)
C
      INCLUDE 'params.inc'
      DIMENSION ERRVEC(NERVEC)
C
C     Internal base vectors.
C     All pointers in the tables are relative to these base vectors, and
C     an offset to the users base vector is maintained to give the
C     user the correct pointer to the users base vector.  They are
C     equivalenced to allow printing as real or integer (see mdprnt).
C
      double precision dmyv(1)
      REAL RMYV(1)
      equivalence (dmyv,rmyv)
c
      DIMENSION MYV(1)
      EQUIVALENCE (MYV,RMYV)
c
      PARAMETER (MAXSIZ=2)
      CHARACTER*1 MYCHAR(MAXSIZ), CFDATA
C      COMMON /EXTCLB/ MYCHAR
C
      LOGICAL INIT, CINIT, DEFER, DEFROK, FILL, CFILL
C
C     Users base vector.  This is an argument in the entry point
C     MDINIT.
C
      DIMENSION USERV(1)
      CHARACTER*1 UCHAR(1)
C
      PARAMETER (NCOLP=132)
      CHARACTER*8  NAME1, NAME2, NAMET
      CHARACTER*(*)  TNAME1, TNAME2
      DIMENSION NVOIDS(2), NNAMES(2)
      SAVE
C
C     Data initialization
C
      DATA NVOIDS /2*0/, NNAMES /2*0/
      DATA ERRVEC /NERVEC*0/
      DATA INIT /.FALSE./, CINIT /.FALSE./,
     *   DEFER /.FALSE./, DEFROK /.FALSE./,
     *   FILL /.FALSE./, CFILL /.FALSE./
      DATA EUNIT /0/
      DATA LERROR /0/
      DATA CHRCOL /1/
C
      LDICT = 2 * LPDICT
      LVOID = 2 * LPVOID
C
C     Get the offset between "USERV" and "MYV".  All internal
C     information is in terms of "MYV" since most operations are done
C     with this base address.
C
      MYLOC = IXLNUM (MYV)
      ULOC = IXLNUM (USERV)
      OFFSET = MYLOC - ULOC
C
C     Check to see if deferred mode will be allowed.  This requires
C     a functioning IXLNUM.
C      IF (IXLNUM(NVOIDS(1)) .NE. IXLNUM(NVOIDS(2))) DEFROK = .TRUE.
C ... Using malloc/free -- ignore deferred mode
      DEFROK = .false.
C
C     Was MDINIT previously called?
C
      IF (INIT) THEN
         LASTER = REINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *)'SECOND CALL TO MDINIT NOT ',
     *      'ALLOWED.'
      ELSE
         LASTER = SUCESS
      END IF
      INIT = .TRUE.
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) LERROR = LASTER
      RETURN
C
C***********************************************************************
C
      ENTRY MCINIT (UCHAR)
C
C***********************************************************************
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MCINIT.'
         RETURN
      END IF
C
      MYCLOC = IXLCHR (MYCHAR(1))
      UCLOC = IXLCHR (UCHAR(1))
      TEMP = IXLCHR (MYCHAR(2))
      IF (MYCLOC .NE. TEMP) THEN
C
C        Mixed mode - character and numeric data intermixed.
C
         CALL EXPARM (NAMET, NAMET, I, CHRNUM, I, I)
         COFFST = CHRNUM * MYLOC - UCLOC
      ELSE
C
C        Split mode - character and numeric data separate.
C
         COFFST = 0
         CHRCOL = 2
         CHRNUM = 1
         IF (NNAMES(1) .GT. LPDICT) THEN
            LASTER = DFULL
         ELSE IF (NVOIDS(1) .GT. LPVOID) THEN
            LASTER = VFULL
         ELSE
            LDICT = LPDICT
            LVOID = LPVOID
            LASTER = SUCESS
         END IF
      END IF
C
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) 'ERROR IN MCINIT.'
         END IF
      ELSE
         CINIT = .TRUE.
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCCOMP
C
C***********************************************************************
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT WAS NOT CALLED ',
     *      'BEFORE CALLING MCCOMP.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C***********************************************************************
C
      ENTRY MDCOMP
C
C***********************************************************************
C
C     This section compresses all void space to the end of each block.
C     No attempt is made to shift storage between blocks.
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDCOMP.'
         RETURN
      END IF
C
C ... Using malloc/free -- ignore MCCOMP/MDDOMP Calls

c      CALL MXCOMP (MYV, VOID, LVOID,
c     *   NVOIDS, DPOINT, LDICT, NNAMES, CHRCOL, LASTER)
C
c      IF (LASTER .EQ. SUCESS .AND. CHRCOL .EQ. 2) THEN
c         CALL MYCOMP (MYCHAR, VOID, LVOID,
c     *      NVOIDS, DPOINT, LDICT, NNAMES, CHRCOL, LASTER)
c      END IF
C
      laster = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) 'ERROR IN MDCOMP/MCCOMP.'
         END IF
      END IF
C
      RETURN
C
C***********************************************************************
C
      ENTRY MCDEBG (UNIT)
      ENTRY MDDEBG (UNIT)
C
C***********************************************************************
C
      IF (UNIT .GE. 1) THEN
         EUNIT = UNIT
      ELSE
         EUNIT = 0
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MDDEL (TNAME1)
C
C***********************************************************************
C
C     This entry point deletes an entry from the dictionary and
C     inserts an entry in the void table.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDDEL.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT WAS NOT CALLED ',
     *      'BEFORE CALLING MDDEL.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      CALL MXDEL (NAME1, DICT, DPOINT, LDICT, NNAMES, VOID,
     *   LVOID, NVOIDS, CHRCOL, LASTER, MYLOC)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MDDEL: ', NAME1
      end if
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MDDEL ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCDEL (TNAME1)
C
C***********************************************************************
C
C     This entry point deletes an entry from the dictionary and
C     inserts an entry in the void table.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCDEL.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT WAS NOT CALLED ',
     *      'BEFORE CALLING MCDEL.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      CALL MYDEL (NAME1, DICT, DPOINT, LDICT, NNAMES, VOID,
     *   LVOID, NVOIDS, CHRCOL, LASTER, MYLOC, MYCLOC)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MCDEL: ', NAME1
      end if
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MCDEL ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCEFIX (ENUM, COUNT)
      ENTRY MDEFIX (ENUM, COUNT)
C
C***********************************************************************
C
C     This resets the number of occurances of a particular error.
C     If the requested error number (enum) is out of range, an error
C     is flagged.
C
      IF (ENUM .GE. 1 .AND. ENUM .LE. NERVEC) THEN
         ERRVEC(ENUM) = COUNT
         LASTER = SUCESS
      ELSE
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'ILLEGAL ERROR FLAG NUMBER',
     *      ' SPECIFIED IN MDEFIX.'
         LASTER = BADCOD
      END IF
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCEROR (UNIT)
      ENTRY MDEROR (UNIT)
C
C***********************************************************************
C
C     This section prints the error codes.
C
      CALL MXEROR (UNIT, LASTER, ERRVEC)
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCERPT (ENUM, COUNT)
      ENTRY MDERPT (ENUM, COUNT)
C
C***********************************************************************
C
C     This reports the number of occurances of a particular error.
C     If the requested error number (ENUM) is out of range, a -1
C     is returned.
C
      IF (ENUM .GE. 1 .AND. ENUM .LE. NERVEC) THEN
         COUNT = ERRVEC(ENUM)
      ELSE
         COUNT = -1
      END IF
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCEXEC ()
C
C***********************************************************************
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCEXEC.'
         RETURN
      END IF
C
C***********************************************************************
C
      ENTRY MDEXEC ()
C
C***********************************************************************
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDEXEC.'
         RETURN
      END IF
C
C ... Using malloc/free -- ignore MCEXEC/MDEXEC calls.
      
c      CALL MXEXEC (MYV, MYCHAR, MYLOC, MYCLOC, UCLOC, COFFST,
c     *   OFFSET, DPOINT, LDICT, NNAMES,
c     *   VOID, LVOID, NVOIDS, FILL, FDATA, CFILL, CFDATA, CHRNUM,
c     *   CHRCOL, LASTER)
c      DEFER = .FALSE.
      laster = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MDDOFF ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MDFILL (TFILL)
C
C***********************************************************************
C
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDFILL.'
         RETURN
      END IF
C
      FDATA = TFILL
      FILL = .TRUE.
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCFILL (TNAME1)
C
C***********************************************************************
C
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCFILL.'
         RETURN
      END IF
C
      IF (LEN(TNAME1) .LT. 1) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BAD CHARACTER STRING '//
     *      'INPUT TO MCFILL.'
         RETURN
      END IF
      CFDATA = TNAME1(1:1)
      CFILL = .TRUE.
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MDFIND (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section finds a name in the dictionary and returns the
C     location and length of that vector.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDFIND.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDFIND.'
         RETURN
      END IF
C
C     Find the name in the dictionary.
C
      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
C
      IF (LASTER .EQ. SUCESS) THEN
C
C        Entry was found in dictionary.
C
         IF (DPOINT(ROW,1,2) .LT. 0) THEN
            DPOINT(ROW,1,1) = IXLNUM(NEWLOC)
            NEWLOC = - ULOC
         ELSE
            NEWLOC = DPOINT(ROW,1,1) + OFFSET
         END IF
         NEWLEN = ABS(DPOINT(ROW,1,2))
      END IF
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) 'NAME ', NAME1, ' MDFIND ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCFIND (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section finds a name in the dictionary and returns the
C     location and length of that vector.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCFIND.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCFIND.'
         RETURN
      END IF
C
C     Find the name in the dictionary.
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
C
      IF (LASTER .EQ. SUCESS) THEN
C
C        Entry was found in dictionary.
C
         IF (DPOINT(ROW,CHRCOL,2) .LT. 0) THEN
            DPOINT(ROW,CHRCOL,1) = IXLNUM(NEWLOC)
            NEWLOC = - UCLOC
         ELSE IF (DPOINT(ROW,CHRCOL,3) .EQ. 0) THEN
            NEWLOC = 1
         ELSE
            NEWLOC = (DPOINT(ROW,1,1) - 1) * CHRNUM + 1 + COFFST
         END IF
         NEWLEN = DPOINT(ROW,1,3)
      END IF
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) 'NAME ', NAME1, ' MCFIND ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MDFOFF ()
C
C***********************************************************************
C
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDFOFF.'
         RETURN
      END IF
C
      FILL = .FALSE.
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCFOFF ()
C
C***********************************************************************
C
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCFOFF.'
         RETURN
      END IF
C
      CFILL = .FALSE.
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MDGET (MNGET)
C
C***********************************************************************
C
C     This section processes a request for a contiguous
C     chunk of memory.
C     If a void of sufficient size is available, no action is taken,
C     otherwise a call to the system for the memory is made.
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDGET.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C ... Using malloc/free -- ignore mdget calls.
      
c      CALL MXGET (MYLOC, MNGET, VOID, LVOID, NVOIDS,
c     *   CHRCOL, LASTER, VROW)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MDGET: ', 0
      end if
      RETURN
C
C***********************************************************************
C
      ENTRY MCGET (MNGET)
C
C***********************************************************************
C
C     This section processes a request for a contiguous
C     chunk of memory.
C     If a void of sufficient size is available, no action is taken,
C     otherwise a call to the system for the memory is made.
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCGET.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C ... Using malloc/free -- ignore mdget calls.
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MCGET: ', 0
      end if
      RETURN
C
C***********************************************************************
C
      ENTRY MCGIVE
C
C***********************************************************************
C
C     This section returns unused space at the end of all storage
C     blocks to the system.  It does not return character memory
C     in split mode.
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCGIVE.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C***********************************************************************
C
      ENTRY MDGIVE
C
C***********************************************************************
C
C     This section returns unused space at the end of all storage
C     blocks to the system.
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDGIVE.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C ... Using malloc/free -- ignore mcgive/mdgive calls.
C
c      CALL MXGIVE (MYLOC, DPOINT, LDICT, NNAMES, VOID, LVOID,
c     *   NVOIDS, CHRCOL, LASTER)
C
      RETURN
C
C***********************************************************************
C
      ENTRY MCLAST (ENUM)
      ENTRY MDLAST (ENUM)
C
C***********************************************************************
C
      ENUM = LERROR
      RETURN
C
C***********************************************************************
C
      ENTRY MCLIST (UNIT)
C
C***********************************************************************
C
C     This section lists the internal tables of the data manager.
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCLIST.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C***********************************************************************
C
      ENTRY MDLIST (UNIT)
C
C***********************************************************************
C
C     This section lists the internal tables of the data manager.
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDLIST.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      CALL MXLIST (UNIT, OFFSET, DICT, DPOINT, LDICT, NNAMES,
     *   VOID, LVOID, NVOIDS, CHRCOL)
C
      if (eunit .gt. 0) then
         CALL MXLIST (EUNIT, OFFSET, DICT, DPOINT, LDICT, NNAMES,
     *        VOID, LVOID, NVOIDS, CHRCOL)
      end if
      
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MDLONG (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section expands or shortens a vector.
C     Any nonnegative value of newlen is permissible.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDLONG.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDLONG.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Is the new length sensible?
C
      IF (NEWLEN .LT. 0) THEN
         LASTER = BADLEN
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'NEGATIVE LENGTH REQUEST ',
     *      '(', NEWLEN, ') FOR NAME ', NAME1, ' IN MDLONG.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C  Ensure that the amount of space called for is even --- this
C  guarantees that the memory will align on an eight byte boundary.
C  (this is only a quick fix for the HP750 --- jrr 10/25/91)
C
      dummy = mod( newlen, 2 )
      if ( dummy .ne. 0 ) then
         junk = newlen + ( 2 - dummy )
      else
         junk = newlen
      endif
C
#ifdef SDP
      INUM = junk
#else
      INUM = newlen
#endif
      CALL MXLONG (NAME1, INUM, NEWLOC, MYV, MYCHAR, MYLOC,
     *   MYCLOC, UCLOC, COFFST, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, VOID, LVOID, NVOIDS,
     *   FILL, FDATA, CFILL, CFDATA, CHRNUM, CHRCOL, LASTER)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MDLONG: ', NAME1, INUM, newloc
      end if
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MDLONG ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCLONG (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section expands or shortens a vector.
C     Any nonnegative value of newlen is permissible.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCLONG.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCLONG.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Is the new length sensible?
C
      IF (NEWLEN .LT. 0) THEN
         LASTER = BADLEN
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'NEGATIVE LENGTH REQUEST ',
     *      '(', NEWLEN, ') FOR NAME ', NAME1, ' IN MCLONG.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C  Ensure that the amount of space called for is a multiple of eight
c  ---this guarantees that the memory will align on an eight byte boundary.
C  (this is only a quick fix for the HP750 --- jrr 10/25/91)
C
      dummy = mod( newlen, 8 )
      if ( dummy .ne. 0 ) then
         junk = newlen + ( 8 - dummy )
      else
         junk = newlen
      endif
#ifdef SDP
      INUM = junk
#else
      INUM = newlen
#endif
C
      CALL MYLONG (NAME1, INUM, NEWLOC, MYV, MYCHAR, MYLOC,
     *   MYCLOC, UCLOC, COFFST, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, VOID, LVOID, NVOIDS,
     *   FILL, FDATA, CFILL, CFDATA, CHRNUM, CHRCOL, MAXSIZ, LASTER)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MCLONG: ', NAME1, INUM, newloc
      end if

      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MCLONG ERROR.'
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MDMEMS (NSUA, NSUD, NSUV, NSULV)
C
C***********************************************************************
C
C     This section returns:
C     NSUA  Number of numeric storage units allocated
C     NSUD  Number of numeric storage units deferred
C     NSUV  Number of numeric storage units in voids
C     NSULV Number of numeric storage units in largest void
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDMEMS.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      NSUA = 0
      NSUD = 0
      NSUV = 0
      NSULV= 0
C
      DO 100 I = 1, NNAMES(1)
         IF (DPOINT(I,1,3) .EQ. -1) THEN
            NSUA = NSUA + MAX(DPOINT(I,1,2),0)
            NSUD = NSUD + MIN(DPOINT(I,1,2),0)
         END IF
  100 CONTINUE
C
      DO 110 I = 1, NVOIDS(1)
         NSUV = NSUV + VOID(I,1,2)
         NSULV = MAX ( NSULV, VOID(I,1,2) )
  110 CONTINUE
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MCMEMS (NSUA, NSUD, NSUV, NSULV)
C
C***********************************************************************
C
C     This section returns:
C     NSUA  Number of character storage units allocated
C     NSUD  Number of character storage units deferred
C     NSUV  Number of character storage units in voids
C     NSULV Number of character storage units in largest void
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MCMEMS.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      NSUA = 0
      NSUD = 0
      NSUV = 0
      NSULV= 0
C
      DO 120 I = 1, NNAMES(CHRCOL)
         IF (DPOINT(I,CHRCOL,3) .GE. 0) THEN
            IF (DPOINT(I,CHRCOL,2) .GE. 0) THEN
               NSUA = NSUA + DPOINT(I,CHRCOL,3)
            ELSE
               NSUD = NSUD + DPOINT(I,CHRCOL,3)
            END IF
         END IF
  120 CONTINUE
C
      DO 130 I = 1, NVOIDS(CHRCOL)
         NSUV = NSUV + VOID(I,CHRCOL,2)
         NSULV = MAX ( NSULV, VOID(I,CHRCOL,2) )
  130 CONTINUE
      NSUV = NSUV * CHRNUM
      NSULV = NSULV * CHRNUM
      LASTER = SUCESS
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      RETURN
C
C***********************************************************************
C
      ENTRY MDNAME (TNAME1, TNAME2)
C
C***********************************************************************
C
C     This section renames an array from NAME1 to NAME2.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDNAME.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
      NAMET = TNAME2
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDNAME.'
         RETURN
      END IF
      NAME2 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDNAME.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Find NAME1 in dictionary.
C
      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .EQ. SUCESS) THEN
C
         LOC = DPOINT(ROW,1,1)
         LEN1 = DPOINT(ROW,1,2)
C
C        Delete dictionary entry.
C
         CALL SHFTC (DICT, CHRCOL*LDICT, ROW+1, NNAMES(1), 1)
         CALL SHFTI (DPOINT, LDICT*CHRCOL, 3, ROW+1, NNAMES(1), 1)
         NNAMES(1) = NNAMES(1) - 1
C
C        insert NAME2 with old pointer and length.
C
         CALL MXNSRT (NAME2, LOC, LEN1, DICT, DPOINT, LDICT,
     *      NNAMES, CHRCOL, LASTER)
         IF (LASTER .NE. SUCESS) THEN
            CALL MXNSRT (NAME1, LOC, LEN1, DICT, DPOINT, LDICT,
     *         NNAMES, CHRCOL, TMP)
         END IF
C
      END IF
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *)
     *         ' MDNAME ERROR.', ' FROM: ', NAME1, ' TO: ', NAME2
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MCNAME (TNAME1, TNAME2)
C
C***********************************************************************
C
C     This section renames an array from NAME1 to NAME2.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCNAME.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
      NAMET = TNAME2
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCNAME.'
         RETURN
      END IF
      NAME2 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCNAME.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Find NAME1 in dictionary.
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .EQ. SUCESS) THEN
C
         LOC = DPOINT(ROW,CHRCOL,1)
         LEN1 = DPOINT(ROW,CHRCOL,2)
         LEN2 = DPOINT(ROW,CHRCOL,3)
C
C        Delete dictionary entry.
C
         CALL SHFTC (DICT(1,CHRCOL), CHRCOL*LDICT, ROW+1,
     *      NNAMES(CHRCOL), 1)
         CALL SHFTI (DPOINT(1,CHRCOL,1), LDICT*CHRCOL, 3, ROW+1,
     *      NNAMES(CHRCOL), 1)
         NNAMES(CHRCOL) = NNAMES(CHRCOL) - 1
C
C        insert NAME2 with old pointer and length.
C
         CALL MYNSRT (NAME2, LOC, LEN1, LEN2, DICT, DPOINT, LDICT,
     *      NNAMES, CHRCOL, LASTER)
         IF (LASTER .NE. SUCESS) THEN
            CALL MYNSRT (NAME1, LOC, LEN1, LEN2, DICT, DPOINT, LDICT,
     *         NNAMES, CHRCOL, TMP)
         END IF
C
      END IF
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *)
     *         ' MCNAME ERROR.', ' FROM: ', NAME1, ' TO: ', NAME2
         END IF
      END IF
      RETURN
C
C***********************************************************************
C
      ENTRY MDPRNT (TNAME1, UNIT, TNAME2)
C
C***********************************************************************
C
C     This section prints a named vector as integer or real, as
C     indicated by NAME2.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDPRNT.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
      NAMET = TNAME2
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK DATA TYPE IN MDPRNT.'
         RETURN
      END IF
      NAME2 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDPRNT.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      CALL MXPRNT (NAME1, UNIT, NAME2, MYV, RMYV, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, CHRCOL, NCOLP, LASTER)
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MDPRNT ERROR.'
         END IF
      END IF
C
      RETURN
C
C***********************************************************************
C
      ENTRY MCPRNT (TNAME1, UNIT, NEWLEN)
C
C***********************************************************************
C
C     This section prints a named vector as character.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCPRNT.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCPRNT.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      CALL MYPRNT (NAME1, UNIT, MYCHAR, COFFST, COFFST+UCLOC-MYCLOC,
     *   DICT, DPOINT, LDICT, NNAMES, CHRNUM,
     *   CHRCOL, NCOLP, NEWLEN, LASTER)
C
      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MCPRNT ERROR.'
         END IF
      END IF
C
      RETURN
C
C***********************************************************************
C
      ENTRY MDRSRV (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section reserves space for an array.  If a void is not
C     available, a system call will be made for the space.  Any
C     nonnegative request is permissible.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MDRSRV.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDRSRV.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Check for sensible length.
C
      IF (NEWLEN .LT. 0) THEN
         LASTER = BADLEN
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'NEGATIVE LENGTH REQUEST ',
     *      '(', NEWLEN, ') FOR NAME ', NAME1, ' IN MDRSRV.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C  Ensure that the amount of space called for is even --- this
C  guarantees that the memory will align on an eight byte boundary.
C  (this is only a quick fix for the HP750 --- jrr 10/25/91)
C
      dummy = mod( newlen, 2 )
      if ( dummy .ne. 0 ) then
         junk = newlen + ( 2 - dummy )
      else
         junk = newlen
      endif
#ifdef SDP
      INUM = junk
#else
      INUM = newlen
#endif
C
      CALL MXRSRV (MYV, NAME1, INUM, NEWLOC, MYLOC, OFFSET,
     *   VOID, LVOID, NVOIDS,
     *   DICT, DPOINT, LDICT, NNAMES, CHRCOL, DEFER, FILL, FDATA,
     *   LASTER)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MDRSRV: ', NAME1, INUM, newloc
      end if

      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MDRSRV ERROR.'
         END IF
      END IF
C
      RETURN
C
C***********************************************************************
C
      ENTRY MCRSRV (TNAME1, NEWLOC, NEWLEN)
C
C***********************************************************************
C
C     This section reserves space for an array.  If a void is not
C     available, a system call will be made for the space.  Any
C     nonnegative request is permissible.
C
      NAMET = TNAME1
      CALL EXUPCS (NAMET)
      CALL STRIPB (NAMET, LEFT, RIGHT)
      IF (LEFT .GT. RIGHT) THEN
         LASTER = BADNAM
         LERROR = LASTER
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'BLANK NAME IN MCRSRV.'
         RETURN
      END IF
      NAME1 = NAMET(LEFT:RIGHT)
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCRSRV.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C     Check for sensible length.
C
      IF (NEWLEN .LT. 0) THEN
         LASTER = BADLEN
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'NEGATIVE LENGTH REQUEST ',
     *      '(', NEWLEN, ') FOR NAME ', NAME1, ' IN MCRSRV.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
C  Ensure that the amount of space called for is a multiple of eight
c  ---this guarantees that the memory will align on an eight byte boundary.
C  (this is only a quick fix for the HP750 --- jrr 10/25/91)
C
      dummy = mod( newlen, 8 )
      if ( dummy .ne. 0 ) then
         junk = newlen + ( 8 - dummy )
      else
         junk = newlen
      endif
#ifdef SDP
      INUM = junk
#else
      INUM = newlen
#endif
C
      CALL MYRSRV (MYCHAR, NAME1, INUM, NEWLOC, MYLOC, MYCLOC,
     *   UCLOC, OFFSET, COFFST,
     *   VOID, LVOID, NVOIDS,
     *   DICT, DPOINT, LDICT, NNAMES, CHRCOL, CHRNUM, DEFER,
     *   CFILL, CFDATA, MAXSIZ,
     *   LASTER)
C
      if (eunit .gt. 0) then
         write (eunit,*) 'MCRSRV: ', NAME1, INUM, newloc
      end if

      ERRVEC(LASTER) = ERRVEC(LASTER) + 1
      IF (LASTER .NE. SUCESS) THEN
         LERROR = LASTER
         IF (EUNIT .GT. 0) THEN
            CALL MXEROR (EUNIT, LASTER, ERRVEC)
            WRITE (EUNIT, *) ' MCRSRV ERROR.'
         END IF
      END IF
C
      RETURN
C
C***********************************************************************
C
      ENTRY MDSTAT (MNERRS, MNUSED)
C
C***********************************************************************
C
C     This section returns the total number of errors and the total
C     amount of numeric storage assigned to arrays.
C
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      MNERRS = 0
      DO 140 I = 2, NERVEC
         MNERRS = MNERRS + ERRVEC(I)
  140 CONTINUE
      MNUSED = 0
      DO 150 I = 1, NNAMES(1)
         IF (DPOINT(I,1,3) .EQ. -1) THEN
            MNUSED = MNUSED + ABS(DPOINT(I,1,2))
         END IF
  150 CONTINUE
      if (eunit .gt. 0) then
         write (eunit,*) 'MDSTAT: ',MNERRS, MNUSED
      end if
      RETURN
C
C***********************************************************************
C
      ENTRY MCSTAT (MNERRS, MNUSED)
C
C***********************************************************************
C
C     This section returns the total number of errors and the total
C     amount of character storage assigned to arrays.
C
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      MNERRS = 0
      DO 160 I = 2, NERVEC
         MNERRS = MNERRS + ERRVEC(I)
  160 CONTINUE
      MNUSED = 0
      DO 170 I = 1, NNAMES(CHRCOL)
         MNUSED = MNUSED + MAX (0, DPOINT(I,CHRCOL,3))
  170 CONTINUE
      if (eunit .gt. 0) then
         write (eunit, *) 'MCSTAT: ', MNERRS, MNUSED
      end if
      RETURN
C
C***********************************************************************
C
      ENTRY MCWAIT ()
C
C     Was MCINIT previously called?
C
      IF (.NOT. CINIT) THEN
         LASTER = NOINIT
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MCINIT NOT CALLED BEFORE ',
     *      'CALLING MCWAIT.'
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         LERROR = LASTER
         RETURN
      END IF
C
      ENTRY MDWAIT ()
C
C***********************************************************************
C
C
C     Was MDINIT previously called?
C
      IF (.NOT. INIT) THEN
         LASTER = NOINIT
         ERRVEC(LASTER) = ERRVEC(LASTER) + 1
         IF (EUNIT .GT. 0) WRITE (EUNIT, *) 'MDINIT NOT CALLED BEFORE ',
     *      'CALLING MDWAIT.'
         RETURN
      END IF
C
      IF (DEFROK) DEFER = .TRUE.
      ERRVEC(SUCESS) = ERRVEC(SUCESS) + 1
      RETURN
C
      END
