C     Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
C     of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C     NTESS, the U.S. Government retains certain rights in this software.
C     
C     See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PREB1 (OPTION, NOUT, IEL0, NLISEL, LISEL,
     &     NLINK, NATR, LINK, ATRIB, natrdm,
     *     MAPN, DOMAPN, MAPE, DOMAPE)
C=======================================================================

C     --*** PREB1 *** (BLOT) Display database element block information
C     --
C     --PREB1 displays the element block internal information (either the
C     --connectivity or the attributes.
C     --
C     --Parameters:
C     --   OPTION - IN - '*' to print all, else print options:
C     --      'C' to print connectivity
C     --      'A' to print attributes
C     --   NOUT - IN - the output file, <=0 for standard
C     --   IEL0 - IN - the element offset of the elements in the block
C     --   NLISEL - IN - the number of selected elements
C     --   LISEL - IN - the indices of the selected elements
C     --   NLINK - IN - the number of nodes per element
C     --   NATR - IN - the number of attributes per element
C     --   LINK - IN - the connectivity array for this block
C     --   ATRIB - IN - the attribute array for this block

      include 'exp_dbnums.blk'
      CHARACTER*(*) OPTION
      INTEGER LISEL(*)
      INTEGER LINK(NLINK,*)
      REAL ATRIB(natrdm,*)
      INTEGER MAPN(*), MAPE(*)
      CHARACTER*32 STRA
      LOGICAL DOMAPN, DOMAPE

      LOGICAL ALLSAM
      LOGICAL DOCONN, DOATR

      IF (NLISEL .LE. 0) RETURN

      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
     &     .AND. (NLINK .GT. 0)
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))
     &     .AND. (NATR .GT. 0)
      IF (.NOT. (DOCONN .OR. DOATR)) RETURN

      if (doconn) then
         if (domapn) then
            stra = '(Global Node IDs)'
         else
            stra = '(Local Node IDs)'
         end if

         if (nout .gt. 0) then
            if (domape) then
               WRITE (NOUT, 10000) 'Connectivity', STRA(:LENSTR(STRA))
            else
               WRITE (NOUT, 10001) 'Connectivity', STRA(:LENSTR(STRA))
            end if
         else
            if (domape) then
               WRITE (*, 10000) 'Connectivity', STRA(:LENSTR(STRA))
            else
               WRITE (*, 10001) 'Connectivity', STRA(:LENSTR(STRA))
            end if
         end if
         nprnt = 0
         do ix=1, numel - iel0
            IEL = LISEL(IX)
            if (iel .eq. 0) cycle
            if (domape) then
               ID = MAPE(IEL)
            else
               id = iel
            end if
            NE = IEL - IEL0
            IF (NOUT .GT. 0) THEN
               if (domapn) then
                  WRITE (NOUT, 10010, IOSTAT=IDUM)
     &                 NE, ID, (MAPN(LINK(I,NE)), I=1,NLINK)
               else
                  WRITE (NOUT, 10010, IOSTAT=IDUM)
     &                 NE, ID, (LINK(I,NE), I=1,NLINK)
               end if
            ELSE
               if (domapn) then
                  WRITE (*, 10010, IOSTAT=IDUM)
     &                 NE, ID, (MAPN(LINK(I,NE)), I=1,NLINK)
               else
                  WRITE (*, 10010, IOSTAT=IDUM)
     &                 NE, ID, (LINK(I,NE), I=1,NLINK)
               end if
            END IF
            nprnt = nprnt + 1
            if (nprnt .ge. nlisel) exit
         end do
      END IF

      if (doatr) then
C     ... See if all attributes are the same
         allsam = .true.
         ifrst = 0
         do i=1, numel
            if (lisel(i) .gt. 0) then
               ifrst = i
               exit
            end if
         end do

         if (ifrst .eq. 0) return

         do ia = 1, natr
            iel = lisel(ifrst)-iel0
            aval = atrib(ia, iel)
            do ix = ifrst+1, numel - iel0
               iel = lisel(ix)
               if (iel .eq. 0) cycle
               ne = iel - iel0
               if (aval .ne. atrib(ia, ne)) then
                  allsam = .false.
                  go to 190
               end if
            end do
         end do

C     ... Print either all values (if not all same), or the common values
 190     continue
         if (allsam) then
            ne = lisel(ifrst) - iel0
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10025, IOSTAT=IDUM)
     &              (ATRIB(I,NE), I=1,NATR)
            ELSE
               WRITE (*, 10025, IOSTAT=IDUM)
     &              (ATRIB(I,NE), I=1,NATR)
            END IF
         else
            if (nout .gt. 0) then
               WRITE (NOUT, 10000) 'Attributes'
            else
               WRITE (*, 10000) 'Attributes'
            end if
            DO IX = ifrst, numel - iel0
               IEL = LISEL(IX)
               if (iel .eq. 0) cycle
               NE = IEL - IEL0
               IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10020, IOSTAT=IDUM)
     &                 NE, IEL, (ATRIB(I,NE), I=1,NATR)
               ELSE
                  WRITE (*, 10020, IOSTAT=IDUM)
     &                 NE, IEL, (ATRIB(I,NE), I=1,NATR)
               END IF
            end do
         end if
      end if
      
      RETURN
      
10000 FORMAT (1X, '  Blk Index  Global ID       ', A, '  ', A)
10001 FORMAT (1X, '  Blk Index   Local ID       ', A, '  ', A)
10010 FORMAT (1X, I11, I11, 5X, 8I11, :, /,
     &     (18X, 8I11))
10020 FORMAT (1X, I11, I11, 3X, 4 (1X, 1pE11.4), :, /,
     &     (16X, 4 (2X, 1pE11.4)))
10025 FORMAT (2X, 'All attributes are equal. Values are: ',
     &     (3 (1x, 1pE11.4)), :, /,
     &     (3X, 6 (1X, 1pE11.4)))
      END
