c @HEADER
c *****************************************************************************
c  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
c
c Copyright 2012 NTESS and the Zoltan contributors.
c SPDX-License-Identifier: BSD-3-Clause
c *****************************************************************************
c @HEADER
      subroutine mmread(iunit,rep,field,symm,rows,cols,nnz,nnzmax,
     *                 indx,jndx,ival,rval,cval)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c This routine will read data from a matrix market formatted file.
c The data may be either sparse coordinate format, or dense array format.
c
c The unit iunit must be open, and the file will be rewound on return.
c
c 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
c 18-Oct-96   Change in routine name to match C and Matlab routines.
c 30-Oct-96   Bug fixes in mmio.f:
c                  -looping for comment lines
c                  -fixed non-ansi zero stringlength
c                  -incorrect size calculation for skew-symmetric arrays
c 	      Other changes in mmio.f:
c                  -added integer value parameter to calling sequences  
c                  -enforced proper count in size info line
c                  -added routine to count words in string (countwd)
c            (Thanks to G.P.Leendetse and H.Oudshoom for their review
c             of the initial version and suggested fixes.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c   Arguments:
c
c   name     type      in/out description
c   ---------------------------------------------------------------
c         
c   iunit    integer     in   Unit identifier for the file
c                             containing the data to be read.
c                             Must be open prior to call.
c                             Will be rewound on return.
c         
c   rep     character*10 out  Matrix Market 'representation' 
c                             indicator. On return:
c                      
c                                coordinate   (for sparse data)
c                                array        (for dense data)
c                                elemental    (to be added)    
c                                   
c   field   character*7  out  Matrix Market 'field'. On return:
c                                   
c                                real 
c                                complex
c                                integer
c                                pattern
c                                   
c   symm    character*19 out  Matrix Market 'field'. On return:
c                                   
c                                symmetric
c                                hermitian
c                                skew-symmetric
c                                general          
c         
c   rows     integer     out  Number of rows in matrix.
c         
c   cols     integer     out  Number of columns in matrix.
c         
c   nnz      integer     out  Number of nonzero entries required to
c                             store matrix.
c         
c   nnzmax   integer     in   Maximum dimension of data arrays.
c         
c   indx     integer(nnz)out  Row indices for coordinate format.
c                             Undefined for array format.
c         
c   jndx     integer(nnz)out  Column indices for coordinate format.
c                             Undefined for array format.
c         
c   ival     integer(nnz) out Integer data (if applicable, see 'field')
c         
c   rval     double(nnz) out  Real data (if applicable, see 'field')
c         
c   cval     complex(nnz)out  Complex data (if applicable, see 'field')
c         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c Declarations:
c
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      double precision rpart,ipart
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, nnzmax, iunit
      integer count
      character mmhead*15
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
c
c Read header line and check validity:
c
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
c
c Convert type code to lower case for easier comparisons:
c
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
         stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
c
c Test input qualifiers:
c
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' .and. 
     *    field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 9000
c
c Read through comment lines, ignoring content:
c
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
c
c Just read a non-comment.
c   Now, back up a line, and read for first int, and back up
c   again. This will set pointer to just before apparent size
c   info line.
c   Before continuing with free form input, count the number of
c   words on the size info line to ensure there is the right amount
c   of info (2 words for array matrices, 3 for coordinate matrices).
c
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
c
c   Correct number of words are present, now back up and read them.
c
      backspace (iunit)
c
      if ( rep .eq. 'coordinate' ) then 
c
c Read matrix in sparse coordinate format
c
        read (iunit,fmt=*) rows,cols,nnz
c
c Check to ensure adequate storage is available
c
        if ( nnz .gt. nnzmax ) then
          print *,'insufficent array lengths for matrix of ',nnz,
     *            ' nonzeros.' 
          print *,'resize nnzmax to at least ',nnz,'. (currently ',
     *            nnzmax,')'
          stop
        endif
c
c Read data according to data type (real,integer,complex, or pattern)
c
        if ( field .eq. 'integer' ) then
          do 30 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),ival(i)
 30       continue
        elseif ( field .eq. 'real' ) then
          do 35 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rval(i)
 35       continue
        elseif ( field .eq. 'complex' ) then
          do 40 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 40       continue
        elseif ( field .eq. 'pattern' ) then
          do 50 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i)
 50       continue 
        else 
           print *,'''',field,''' data type not recognized.'
           stop
        endif
        rewind(iunit)
        return
c
      elseif ( rep .eq. 'array' ) then
c
c Read matrix in dense column-oriented array format
c
        read (iunit,fmt=*) rows,cols
c
c Check to ensure adequate storage is available
c
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnzreq = (rows*cols - rows)/2 + rows
          nnz = nnzreq
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnzreq = (rows*cols - rows)/2 
          nnz = nnzreq
        else
          nnzreq = rows*cols
          nnz = nnzreq
        endif
        if ( nnzreq .gt. nnzmax ) then
          print *,'insufficent array length for ',rows, ' by ',
     *             cols,' dense ',symm,' matrix.'
          print *,'resize nnzmax to at least ',nnzreq,'. (currently ',
     *             nnzmax,')'
          stop
        endif
c
c Read data according to data type (real,integer,complex, or pattern)
c
        if ( field .eq. 'integer' ) then
          do 60 i=1,nnzreq
            read (iunit,fmt=*,end=4000) ival(i)
 60      continue
        elseif ( field .eq. 'real' ) then
          do 65 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rval(i)
 65      continue
        elseif ( field .eq. 'complex' ) then
          do 70 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 70      continue
        else
           print *,'''pattern'' data not consistant with type ''array'''
           stop
        endif
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
c
c Various error conditions:
c
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data lines found.'
      stop
 3000 print *,'Size info inconsistant with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistant with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 4000 print *,'Premature end-of-file.'
      print *,'Check that the data file contains ',nnz,
     *        ' lines of  i,j,[val] data.'
      print *,'(it appears there are only ',i,' such lines.)'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c End of subroutine mmread
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      subroutine mminfo(iunit,rep,field,symm,rows,cols,nnz)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c This routine will read header information from a Matrix Market 
c formatted file.  
c
c The unit iunit must be open, and the file will be rewound on return.
c
c 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
c 18-Oct-96   Change in routine name to match C and Matlab routines.
c 30-Oct-96   Bug fixes in mmio.f:
c                  -looping for comment lines
c                  -fixed non-ansi zero stringlength
c                  -incorrect size calculation for skew-symmetric arrays
c 	      Other changes in mmio.f:
c                  -added integer value parameter to calling sequences  
c                  -enforced proper count in size info line
c                  -added routine to count words in string (countwd)
c            (Thanks to G.P.Leendetse and H.Oudshoom for their review
c             of the initial version and suggested fixes.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c   Arguments:
c
c   name     type      in/out description
c   ---------------------------------------------------------------
c        
c   iunit  integer     in   Unit identifier for the open file
c                             containing the data to be read.
c        
c   rep     character*10 out  Matrix Market 'representation' 
c                             indicator. On return:
c                      
c                                coordinate   (for sparse data)
c                                array        (for dense data)
c                                elemental    (to be added)    
c                                   
c   field   character*7  out  Matrix Market 'field'. On return:
c                                   
c                                real 
c                                complex
c                                integer
c                                pattern
c                                   
c   symm    character*19 out  Matrix Market 'field'. On return:
c                                   
c                                symmetric
c                                hermitian
c                                skew-symmetric
c                                general          
c         
c   rows     integer     out  Number of rows in matrix.
c        
c   cols     integer     out  Number of columns in matrix.
c        
c   nnz      integer     out  Number of nonzero entries required to store 
c                             the matrix.
c        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c Declarations:
c
      integer i, rows, cols, nnz, iunit
      integer count
      character mmhead*14
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
c
c Read header line and check validity:
c
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
c
c Parse words from header line:
c
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
c
c Convert type code to upper case for easier comparisons:
c
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
        stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
c
c Test input qualifiers:
c
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' .and. 
     *    field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 9000
c
c Read through comment lines, ignoring content:
c
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
c KDDKDD Changed max number of comment lines j from 2 to 50, as "do 10" loop
c KDDKDD wasn't working with j=2
      j = 50  
      do 10 i=1,j
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        j = j + 1
 10   continue
 20   continue
c
c Just read a non-comment.
c   Now, back up a line, and read for first int, and back up
c   again. This will set pointer to just before apparent size
c   info line.
c   Before continuing with free form input, count the number of
c   words on the size info line to ensure there is the right amount
c   of info (2 words for array matrices, 3 for coordinate matrices).
c
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
c
c   Correct number of words are present, now back up and read them.
c
      backspace (iunit)
c
      if ( rep .eq. 'coordinate' ) then 
c
c Read matrix in sparse coordinate format
c
        read (iunit,fmt=*) rows,cols,nnz
c
c Rewind before returning 
c
        rewind(iunit)
        return
c
      elseif ( rep .eq. 'array' ) then
c
c Read matrix in dense column-oriented array format
c
        read (iunit,fmt=*) rows,cols
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnz = (rows*cols - rows)/2 + rows
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnz = (rows*cols - rows)/2 
        else
          nnz = rows*cols
        endif
c
c Rewind before returning 
c
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
c
c Various error conditions:
c
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data found.'
      stop
 3000 print *,'Size info inconsistant with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistant with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
      stop
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c End of subroutine mmread 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      subroutine mmwrite(ounit,rep,field,symm,rows,cols,nnz,
     *                    indx,jndx,ival,rval,cval)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c This routine will write data to a matrix market formatted file.
c The data may be either sparse coordinate format, or dense array format.
c
c The unit ounit must be open.
c
c 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
c 18-Oct-96   Change in routine name to match C and Matlab routines.
c 30-Oct-96   Bug fixes in mmio.f:
c                  -looping for comment lines
c                  -fixed non-ansi zero stringlength
c                  -incorrect size calculation for skew-symmetric arrays
c 	      Other changes in mmio.f:
c                  -added integer value parameter to calling sequences  
c                  -enforced proper count in size info line
c                  -added routine to count words in string (countwd)
c            (Thanks to G.P.Leendetse and H.Oudshoom for their review
c             of the initial version and suggested fixes.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c   Arguments:
c
c   name     type      in/out description
c   ---------------------------------------------------------------
c         
c   ounit  integer     in   Unit identifier for the file
c                             to which the data will be written.
c                             Must be open prior to call.
c         
c   rep     character*   in   Matrix Market 'representation' 
c                             indicator. Valid inputs:
c                      
c                                coordinate   (for sparse data)
c                                array        (for dense data)
c                               *elemental*    (to be added)    
c                                   
c   field   character*   in   Matrix Market 'field'. Valid inputs:
c                                   
c                                real 
c                                complex
c                                integer
c                                pattern (not valid for dense arrays)
c                                   
c   symm    character*   in   Matrix Market 'field'. Valid inputs:
c                                   
c                                symmetric
c                                hermitian
c                                skew-symmetric
c                                general          
c         
c   rows     integer     in   Number of rows in matrix.
c         
c   cols     integer     in   Number of columns in matrix.
c         
c   nnz      integer     in   Number of nonzero entries in matrix.
c                             (rows*cols for array matrices)
c         
c   indx     integer(nnz)in   Row indices for coordinate format.
c                             Undefined for array format.
c         
c   jndx     integer(nnz)in   Column indices for coordinate format.
c                             Undefined for array format.
c         
c   ival     integer(nnz) in  Integer data (if applicable, see 'field')
c         
c   rval     double(nnz) in   Real data (if applicable, see 'field')
c         
c   cval     complex(nnz)in   Complex data (if applicable, see 'field')
c         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c Declarations:
c
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, ounit
      character*(*)rep,field,symm
c
c Test input qualifiers:
c
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' )
     *   go to 1000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' .and. 
     *    field .ne. 'pattern') go to 2000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. 
     *    field .ne. 'real' .and. field .ne. 'complex' ) go to 3000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and.
     *    symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric')
     *   go to 4000
c
c Write header line:
c
      write(unit=ounit,fmt=5)rep,' ',field,' ',symm
 5    format('%%MatrixMarket matrix ',11A,1A,8A,1A,20A)
c
c Write size information:
c
      if ( rep .eq. 'coordinate' ) then
         nnzreq=nnz
         write(unit=ounit,fmt=*) rows,cols,nnz
         if ( field .eq. 'integer' ) then
            do 10 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),ival(i)
 10         continue
         elseif ( field .eq. 'real' ) then
            do 20 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),rval(i)
 20         continue
         elseif ( field .eq. 'complex' ) then
            do 30 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),
     *                                real(cval(i)),aimag(cval(i))
 30         continue
         else
c        field .eq. 'pattern' 
            do 40 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i)
 40         continue
         endif
      else
c        rep .eq. 'array'
         if ( symm .eq. 'general' ) then
           nnzreq = rows*cols
         elseif ( symm .eq. 'symmetric' .or. 
     *            symm .eq. 'hermitian' ) then
           nnzreq = (rows*cols - rows)/2 + rows
         else 
c        symm .eq. 'skew-symmetric' 
           nnzreq = (rows*cols - rows)/2 
         endif
         write(unit=ounit,fmt=*)rows,cols
         if ( field .eq. 'integer' ) then
            do 50 i=1,nnzreq
               write(unit=ounit,fmt=*)ival(i)
 50         continue
         elseif ( field .eq. 'real' ) then
            do 60 i=1,nnzreq
               write(unit=ounit,fmt=*)rval(i)
 60         continue
         else
c        field .eq. 'complex' 
            do 70 i=1,nnzreq
               write(unit=ounit,fmt=*)real(cval(i)),aimag(cval(i))
 70         continue
         endif
      endif
      return
c
c Various errors
c
 1000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 2000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 3000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 4000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c End of subroutine mmwrite
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      subroutine lowerc(string,pos,len)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c Convert uppercase letters to lowercase letters in string with
c starting postion pos and length len.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer pos, len
      character*(*) string

      character*26 lcase, ucase
      save lcase,ucase
      data lcase/'abcdefghijklmnopqrstuvwxyz'/
      data ucase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do 10 i=pos,len
        k = index(ucase,string(i:i))
        if (k.ne.0) string(i:i) = lcase(k:k)
 10   continue
      return
      end

      subroutine getwd(word,string,slen,start,next,wlen)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Getwd extracts the first  word from string starting
c     at position start.  On return, next is the position
c     of the blank which terminates the word in string.   
c     If the found word is longer than the allocated space
c     for the word in the calling program, the word will be 
c     truncated to fit.
c     Count is set to the length of the word found.
c     
c 30-Oct-96   Bug fix: fixed non-ansi zero stringlength
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer slen, start, next, begin, space, wlen
      character*(*) word
      character*(*) string

      begin = start
      do 5 i=start,slen
         space = index(string(i:slen),' ')
         if ( space .gt. 1) then
            next = i+space-1
            go to 100
         endif
         begin=begin+1
 5    continue
 100  continue
      wlen=next-begin
      if ( wlen .le. 0 ) then
        wlen = 0
        word = ' '
        return
      endif
      word=string(begin:begin+wlen)
      return
      end

      subroutine countwd(string,slen,start,count)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Countwd counts the number of words in string starting
c     at position start.  On return, count is the number of words.
c 30-Oct-96   Routine added
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      character*(*) string
      integer slen, start, next, wordlength, count
      character tmp2*2

      count = 0
      next = 1
 10   call getwd(tmp2,string,1024,next,next,wordlength)
      if ( wordlength .gt. 0 ) then
         count = count + 1
         go to 10
      endif
      return
      end
