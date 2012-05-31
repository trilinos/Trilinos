C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C 
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
C 

C $Id: pltfnt.f,v 1.1 1993/07/16 16:48:09 gdsjaar Exp $ 
C $Log: pltfnt.f,v $
C Revision 1.1  1993/07/16 16:48:09  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTFNT(FILENM)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      CHARACTER*(*) FILENM
      CHARACTER*512 LOCFIL
      character*512 ACCESS
      CHARACTER*132 TLINE
      CHARACTER*10 IERROR
      INTEGER CHRLEN

      CALL CPUNAL(IUNIT)
      LOCFIL = FILENM
      L = CHRLEN(LOCFIL)
      IF (L.EQ.0) THEN
         RETURN

      END IF

      OPEN (UNIT=IUNIT,FILE=LOCFIL(1:L),FORM='unformatted',STATUS='old',
     *  IOSTAT=IOS)
      if (ios .ne. 0) then
C ... See if environment variable ACCESS is set.  If so, then
C     search for file in "$ACCESS/bin/LOCFIL"        
        call getenv('ACCESS', ACCESS)
        la = chrlen(ACCESS)
        if (la .ne. 0) then
          locfil = ACCESS(:la) // '/bin/' // filenm(:l)
          l = chrlen(locfil)
          OPEN (UNIT=IUNIT,FILE=LOCFIL(1:L),FORM='unformatted',
     *      STATUS='old', ERR=10, IOSTAT=IOS)
        else
          go to 10
        end if
      end if
      
      READ (IUNIT) NOVECT,TEXTP(39),TEXTP(40)
      IF (NOVECT.GT.2300) THEN
         CALL PLTFLU
         TLINE = 'Too many vectors in font file "'//LOCFIL(1:L)//
     *           '"; no font read in.'
         CALL SIORPT('PLTFNT',TLINE,2)
         RETURN

      END IF

      READ (IUNIT,ERR=30,IOSTAT=IOS) (IDEX(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (NVECT(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (XSIZE(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (YSIZE(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X0(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y0(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X1(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y1(I,1),I=1,NOVECT)
      READ (IUNIT) NOVECT,TEXTP(39),TEXTP(40)
      IF (NOVECT.GT.2300) THEN
         CALL PLTFLU
         TLINE = 'Too many vectors in font file "'//LOCFIL(1:L)//
     *           '"; no font read in.'
         CALL SIORPT('PLTFNT',TLINE,2)
         RETURN

      END IF

      READ (IUNIT,ERR=30,IOSTAT=IOS) (IDEX(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (NVECT(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (XSIZE(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (YSIZE(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X0(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y0(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X1(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y1(I,2),I=1,NOVECT)
      CLOSE (IUNIT)
      CALL CPUNDE(IUNIT)
      RETURN

   30 CONTINUE
      CALL CHRIC(IOS,IERROR,LI)
      CALL PLTFLU
      TLINE = 'I/O error '//IERROR(1:LI)//' in reading font file "'//
     *        LOCFIL(1:L)//'"; no font read in.'
      CALL SIORPT('PLTFNT',TLINE,2)
      RETURN

   10 CONTINUE
      CALL CHRIC(IOS,IERROR,LI)
      CALL PLTFLU
      TLINE = 'Open error '//IERROR(1:LI)//' in reading font file "'//
     *        LOCFIL(1:L)//'"; no font read in.'
      CALL SIORPT('PLTFNT',TLINE,2)
      RETURN

      END
