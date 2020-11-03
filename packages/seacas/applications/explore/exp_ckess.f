C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CKESS (NUMESS, LESSEL, LESSNL, NUMEL, NUMNP,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &   LTEESS, LTSESS, FACESS, NSCR, ICHECK, RCHECK, NDIM,
     *   MAPEL, MAPND)
C=======================================================================

C   --*** CKESS *** (EXPLORE) Check database element side sets
C   --
C   --CKESS checks the element side set information.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NUMESS - IN - the number of element side sets
C   --   LESSEL - IN - the number of elements for all sets
C   --   LESSNL - IN - the number of nodes for all sets
C   --   NUMEL - IN - the number of elements
C   --   NUMNP - IN - the number of nodes
C   --   IDESS - IN - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTSESS - IN - the element faces for all sets
C   --   FACESS - IN - the distribution factors for all sets
C   --   NSCR - SCRATCH - size = LESSEL (not all used)
C   --   ICHECK - SCRATCH - size = MAX (NUMESS, LESSEL, LESSNL)
C   --   RCHECK - SCRATCH - size = NUMNP

      include 'exodusII.inc'
      INCLUDE 'exp_dbase.blk'

      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL FACESS(*)
      INTEGER ICHECK(*)
      INTEGER NSCR(*)
      REAL RCHECK(*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)
      LOGICAL DIDHEAD, ALLSAM

      CHARACTER*1024 STRA

C   --Check for unique identifier

      DO 100 IESS = 1, NUMESS
         IF (LOCINT (IDESS(IESS), IESS-1, IDESS) .GT. 0) THEN
            CALL INTSTR (1, 0, IDESS(IESS), STRA, LSTRA)
            CALL PRTERR ('WARNING', 'Element side set ID '
     &         // STRA(:LSTRA) // ' is not unique')
         END IF
  100 CONTINUE

C   --Check number of elements in element side sets

      NESS = 0
      DO 110 IESS = 1, NUMESS
         NESS = MAX (NESS, IXEESS(IESS) + NEESS(IESS) - 1)
  110 CONTINUE

      IF (NESS .NE. LESSEL) THEN
         CALL PRTERR ('WARNING', 'Maximum element index'
     &      // ' in all element side sets does not match total')
      END IF

C   --Check all elements in element side sets are within element range

      CALL CHKRNG (LTEESS, LESSEL, NUMEL, NZERO, NERR)
      IF (NERR .GT. 0) THEN
         CALL PRTERR ('FATAL',
     &      'Element side set element ids are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
         CALL PRTERR ('FATAL',
     &      'Element side set element ids are zero')
      END IF

C   --Check all element faces in element side sets are within range
C ... Since we don't know (or don't want to expend the effort...) the
C     the number of faces for each element, we assume that the maximum
C     number of faces is 4 for 2D and 6 for 3D
      CALL CHKRNG (LTSESS, LESSEL, 2*NDIM, NZERO, NERR)
      IF (NERR .GT. 0) THEN
        CALL PRTERR ('FATAL',
     &    'Element side set faces are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
        CALL PRTERR ('FATAL',
     &    'Element side set faces are zero')
      END IF

C ... Check for duplicate element/sides in a sideset. This causes
C     problems with some analysis codes
      do iess = 1, numess
        call iniint(numel, 0, icheck)
        nel = neess(iess)
        indx = ixeess(iess)
        do j = 0, nel-1
          iel = lteess(indx+j)
          ifa = ltsess(indx+j)
          if (btest(icheck(iel), ifa)) then
            write (stra, 10000) iel, ifa, idess(iess)
10000       FORMAT('SIDESET ERROR: The element face pair ',I12,'.',I1,
     $        ' is duplicated in sideset ', I12,'.')
            call sqzstr(stra, lstra)
            CALL PRTERR ('WARNING', STRA(:lstra))
          else
            icheck(iel) = ibset(icheck(iel), ifa)
          end if
        end do
      end do

c ... Check that the distribution factor count matches the number of nodes
C     in the sideset...
      do iess = 1, numess
        call exgsp(ndb, idess(iess), nsess, ndfss, ierr)
        if (nness(iess) .ne. ndfss .and. ndfss .gt. 0) then
           write (stra, 10002) idess(iess), ndfss, nness(iess)
10002      FORMAT('SIDESET ERROR: In sideset ', I12,
     *          ' the number of distribution factors (', I12,
     *          ') does not match the sideset node count (', I12, ')')
           call sqzstr(stra, lstra)
           CALL PRTERR ('WARNING', STRA(:lstra))
        end if
        if (ndfss .gt. 0) then
           call exgssc(ndb, idess(iess), nscr, ierr)
           numnod = 0
           do i = 1, neess(iess)
              numnod = numnod + nscr(i)
           end do
           if (ndfss .ne. numnod) then
              write (stra, 10001) idess(iess), ndfss, numnod
10001         FORMAT('SIDESET ERROR: In sideset ', I12,
     *             ' the number of distribution factors (', I12,
     *             ') does not match the computed sideset node count (',
     *             I12, ')')
              call sqzstr(stra, lstra)
              CALL PRTERR ('WARNING', STRA(:lstra))
           endif
        end if
      end do

c ... Check for discontinuous sideset distribution factors on a sideset.
C     That is, if node 42 on side 15 has a different df value than node 42 on side 11.
C     This is allowed for in exodus, but most users want a c1 continuous field defined.
      do iess = 1, numess
        call inirea(numnp, 0.0, rcheck)
        didhead = .false.
        IF (NNESS(IESS) .GT. 0) THEN
          IS = IXNESS(IESS)
          IE = IS + NNESS(IESS) - 1
C     ... See if all values are the same
          val = facess(is)
          allsam = .TRUE.
          do i=is+1, ie
            if (facess(i) .ne. val) then
              allsam = .FALSE.
              go to 90
            end if
          end do
 90       continue

          if (.not. allsam) then
C ... Get the number of df/nodes per face and the nodes on the face
C ... NOTE: facess is contiguous over all sidesets,
C           nscr is the number of nodes/df per face
C           icheck is the nodes for the faces.
            call exgssn(ndb, idess(iess), nscr, icheck, ierr)
            ISE = IXEESS(IESS)
            IEE = ISE + NEESS(IESS) - 1
            IDS = IS
            idf = 1
            idn = 1
            do i=ise, iee
              NDFPE = nscr(idf)
              idf=idf+1
              do j=1,ndfpe
                node = icheck(j+idn-1)
                if (rcheck(node) .ne. 0.0 .and.
     *            rcheck(node) .ne. facess(j+ids-1)) then
                  iel  = lteess(i)
                  isid = ltsess(i)
                  if (.not. didhead) then
                    WRITE (*, 10040, IOSTAT=IDUM) idess(iess)
                    didhead = .true.
                  end if
                  write (stra, 10050) mapel(iel), isid,
     *              mapnd(node), rcheck(node), facess(j+ids-1)
                  CALL PRTERR ('CMDSPEC', stra(:lenstr(stra)))
                else
                  rcheck(node) = facess(j+ids-1)
                endif
              end do
              ids = ids + ndfpe
              idn = idn + ndfpe
            end do
          end if
        END IF
      end do
10040 FORMAT('SIDESET DF CONTINUITY ERRORS For Sideset ',I12)
10050 FORMAT('Element ',I12,', Side ',I1, ', Node ',I12,
     *  ': Previous Value = ',1PE11.4,'  Current Value = ',1PE11.4)

      RETURN
      END
