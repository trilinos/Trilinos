/**   ------------------------------------------------------------
 *    Copyright 2000-2007 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <cstdlib>
#include <stk_util/parallel/mpi_filebuf.hpp>
#include <assert.h>

enum { buffer_default_length = 4096 };
enum { buffer_putback_length =   16 };

/*--------------------------------------------------------------------*/

mpi_filebuf::mpi_filebuf()
  : std::streambuf(),
    comm( MPI_COMM_NULL ),
    comm_root( -1 ),
    comm_root_fp( NULL ),
    comm_output( 0 ),
    comm_buffer( NULL ),
    comm_buffer_len( buffer_default_length ),
    comm_time(0.0)
{}

mpi_filebuf::~mpi_filebuf()
{
  close();
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::set_buffer_length( const size_t len )
{
  // If already open then abort
  if ( NULL != comm_buffer ) return (mpi_filebuf *) NULL ;

  // Wait and verify upon the attempt to open
  comm_buffer_len = buffer_putback_length < len ? len : buffer_putback_length ;

  return this ;
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::open(
	MPI_Comm       communicator ,
  const int            root_processor ,
  const std::ios_base::openmode file_mode ,
  const char * const   file_name )
{
  const double start_time = MPI_Wtime();

  // If already open then abort
  if ( NULL != comm_buffer ) return (mpi_filebuf *) NULL ;

  const int mode =
    ( std::ios::in  == file_mode ) ? 'r' : (
    ( std::ios::out == file_mode ) ? 'w' : (
    ( std::ios::app == file_mode ) ? 'a' : -1 ) );

  int err ;
  int rank ;
  int local, global ;
  int data[3] ;

  // Broadcast the selected root processor and 'C' file mode

  data[0] = root_processor ;
  data[1] = mode ;
  data[2] = comm_buffer_len ;

  if ( MPI_SUCCESS != ( err = MPI_Bcast(data,3,MPI_INT,0,communicator) ) )
    MPI_Abort( communicator , err );

  // Verify that all processors have the same root, mode, and buffer length:

  local = data[0] != root_processor || data[1] != mode || data[2] != (signed) comm_buffer_len ;

  if ( MPI_SUCCESS != ( err =
       MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,communicator) ) )
    MPI_Abort( communicator , err );

  if ( global ) {
    comm_time += MPI_Wtime() - start_time ;
    return (mpi_filebuf *) NULL ;
  }

  //--------------------------------------------------------------------
  // Root processor and mode are consistent.
  // All processors try to allocate buffers and the
  // root processor tries to open the file.

  if ( MPI_SUCCESS != ( err =  MPI_Comm_rank( communicator , &rank ) ) )
    MPI_Abort( communicator , err );

  char * const tmp_buf = (char *) std::malloc( comm_buffer_len );
  std::FILE *       tmp_fp  = NULL ;

  local = tmp_buf == NULL ; // Failed allocation ?

  if ( root_processor == rank && ! local ) {
    tmp_fp = std::fopen( file_name , ( ( ( mode == 'r' ) ? "r" :
				    ( mode == 'w' ) ? "w" : "a" ) ) );
#ifdef REDSTORM_SETVBUF
    if (tmp_fp) {
      if (std::setvbuf(tmp_fp, NULL, _IOFBF, 32768) != 0) {
	std::fclose(tmp_fp);
	tmp_fp = 0;
      }
    }
#endif
    local = NULL == tmp_fp ;
  }

  if ( MPI_SUCCESS != ( err =
       MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,communicator) ) )
    MPI_Abort( communicator , err );

  if ( global ) {
    if ( NULL != tmp_buf ) std::free(   tmp_buf ); // Deallocate
    if ( NULL != tmp_fp  ) std::fclose( tmp_fp );  // Close the file
    comm_time += MPI_Wtime() - start_time ;
    return (mpi_filebuf *) NULL ;
  }

  //--------------------------------------------------------------------
  // All memory allocated and root processor openned the file
  // Update the internal members accordingly.

  comm         = communicator ;
  comm_root    = root_processor ;
  comm_root_fp = tmp_fp ;
  comm_buffer  = tmp_buf ;
  comm_output  = mode != 'r' ;

  // If output then set up put-buffer

  if ( comm_output ) setp( comm_buffer, comm_buffer + comm_buffer_len );

  comm_time += MPI_Wtime() - start_time ;

  return this ;
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::close()
{
  mpi_filebuf * tmp = NULL ;

  if ( NULL != comm_buffer ) {

    flush(); // Flush the buffers

    if ( NULL != comm_root_fp ) std::fclose( comm_root_fp ); // Close the file

    std::free( comm_buffer ); // Free the buffer

    if ( comm_output ) setp(NULL,NULL);
    else               setg(NULL,NULL,NULL);

    // Reset the members:

    comm         = MPI_COMM_NULL ;
    comm_root    = -1 ;
    comm_root_fp = NULL ;
    comm_output  = 0 ;
    comm_buffer  = NULL ;

    tmp = this ;
  }

  return tmp ;
}

/*--------------------------------------------------------------------*/
/* Underflow, a global call.
   Read more data from the root processor's file and
   broadcast it to all processors.
*/

int mpi_filebuf::underflow()
{
  const double start_time = MPI_Wtime();

  if ( NULL != comm_buffer && ! comm_output &&     // Open for read
       ( gptr() == NULL || gptr() >= egptr() ) ) { // valid get buffer


    // Length of the buffer, consistent on all processors
    // Entire buffer is offset to accomodate putbacks

    const size_t size = comm_buffer_len - buffer_putback_length ;
    char * const buf  = comm_buffer     + buffer_putback_length ;

    int nread ;
    int err ;

    // Root processor reads from the file and broadcasts the result

    if ( NULL != comm_root_fp ) nread = std::fread(buf,1,size,comm_root_fp);

    if ( MPI_SUCCESS != ( err =
	 MPI_Bcast( &nread, 1, MPI_INT, comm_root, comm ) ) )
      MPI_Abort(comm,err);

    // If the read is successfull then update the get buffer pointers:

    if ( 0 < nread ) {

      // Broadcast the read buffer to all processors:

      if ( MPI_SUCCESS != ( err =
	   MPI_Bcast( buf, nread, MPI_BYTE, comm_root, comm ) ) )
	MPI_Abort(comm,err);

      // Set the get buffer:

      setg( comm_buffer, buf, buf + nread );

      // Return the next character from the file:

      comm_time += MPI_Wtime() - start_time ;

      return *buf ;
    }
  }

  // Failed: set the get buffer to NULL and return EOF
  setg(NULL, NULL, NULL);

  comm_time += MPI_Wtime() - start_time ;

  return EOF;
}

/*--------------------------------------------------------------------*/
/* Overflow, a local call.
    Output complete lines of data on the root processor.
    Increase the buffer size on all other processors.
*/

int mpi_filebuf::overflow( int c )
{
  if ( NULL != comm_buffer && comm_output ) { // open for write

    // Determine current offset and length:
    char * cur_buffer = comm_buffer ;
    size_t cur_offset = pptr()  - cur_buffer ;
    size_t cur_length = epptr() - cur_buffer ;

    assert( cur_offset <= cur_length /* detecting abuse by 'ostream' */ );

    if ( NULL != comm_root_fp ) {
      if ( std::fwrite(cur_buffer,1,cur_offset,comm_root_fp) != cur_offset ) {
	return EOF ; // Write failed
      }
      cur_offset = 0 ;
    }
    else if ( cur_length <= cur_offset ) {
      // Not root processor, ran out of buffer space and
      // cannot write so increase the buffer size:
      cur_buffer = (char *) std::realloc( cur_buffer , cur_length *= 2 );
    }

    // If buffer is still good then reset the put-buffer

    if ( NULL != cur_buffer ) {

      comm_buffer = cur_buffer ;

      setp( cur_buffer + cur_offset, cur_buffer + cur_length );

      if ( c != EOF ) {

	sputc(c);
	return c;
      }
      else {
	return 0;
      }
    }
  }
  return EOF ;
}

/*--------------------------------------------------------------------*/
/* Send output buffers to root processor and
   write them to the output file.
*/

mpi_filebuf * mpi_filebuf::flush()
{
  const double start_time = MPI_Wtime();

  int result = -1 ; // Failure return value

  if ( NULL != comm_buffer && comm_output ) { // Open for write

    int err ;

    result = 0 ;

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    unsigned int cur_len = pptr() - cur_buf ;

    // Determine the global lengths

    char * recv_buf  = NULL ;
    int  * recv_len  = NULL ;
    int  * recv_disp = NULL ;

    int nproc = 0 ;


//  if ( NULL != comm_root_fp ) {

//  It should no be neccessary to allocate recv_len on non-root
//  nodes, but the MPI_Gatherv on Janus always accesses recv_len
//  even on non-root processors which causes a segmentaion
//  violation if recv_len is set to NULL.

    if ( MPI_SUCCESS != ( err = MPI_Comm_size(comm,&nproc) ) )
      MPI_Abort( comm , err );

    recv_len = (int*) std::malloc( sizeof(int) * nproc );

    if ( NULL == recv_len ) MPI_Abort( comm , MPI_ERR_UNKNOWN );

    for (int j = 0 ; j < nproc ; ++j )
      recv_len[j] = 0;
//  }

    // Gather buffer lengths on the root processor

    if ( MPI_SUCCESS != ( err =
	 MPI_Gather(&cur_len,1,MPI_INT,recv_len,1,MPI_INT,comm_root,comm)))
      MPI_Abort( comm , err );

    // Root processor must allocate enough buffer space:

    if ( NULL != comm_root_fp ) {

      recv_len[ comm_root ] = 0 ; // Don't send to self

      int i ;

      if ( NULL == ( recv_disp = (int*) std::malloc( sizeof(int) * (nproc + 1) ) ) )
	result = -1 ;

      if ( 0 == result ) { // Allocation succeeded

	recv_disp[0] = 0 ;

	for ( i = 0 ; i < nproc ; ++i )
	  recv_disp[i+1] = recv_disp[i] + recv_len[i] ;

	if ( 0 < recv_disp[nproc] ) {
	  if ( NULL == ( recv_buf = (char*) std::malloc( recv_disp[nproc] ) ) )
	    result = -1 ;
	}
	else {
	  result = 1 ; // No need to gather!
	}

	if ( -1 != result ) {

	  // Write the root processor's buffer

	  if ( 0 < cur_len ) {
	    if ( std::fwrite(cur_buf,1,cur_len,comm_root_fp) != cur_len )
	      result = -1 ; // Write failed

	    cur_len = 0 ; // Wrote this buffer
	  }
	}
      }
      std::fflush( comm_root_fp );
    }

    // Root process broadcasts that all is well with the allocation

    if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
      MPI_Abort( comm , err );

    if ( 0 == result ) { // All-is-well, need to gather and write

      // Gather the buffers to the root processor

      if ( MPI_SUCCESS != ( err =
	   MPI_Gatherv(cur_buf,  cur_len,             MPI_BYTE,
		       recv_buf, recv_len, recv_disp, MPI_BYTE,
		       comm_root, comm ) ) )
	MPI_Abort( comm , err );

       // Output the buffers, beginning with 'comm_root'

      if ( NULL != comm_root_fp ) {

	int i ;

	for ( i = 1 ; i < nproc && 0 == result ; ++i ) {
	  const int j   = ( i + comm_root ) % nproc ;
	  const unsigned int len = recv_len[j] ;

	  if ( 0 < len )
	    if ( std::fwrite(recv_buf+recv_disp[j],1,len,comm_root_fp) != len )
	      result = -1 ; // Write failed
	}

	std::fflush( comm_root_fp );
      }

      // Broadcast that the write succeeded

      if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
	MPI_Abort( comm , err );
    }
    else if ( 1 == result ) {
      // Did not need to gather

      result = 0 ;
    }

    // Reset the output buffer

    setp( comm_buffer , epptr() );

    // Clean up allocated memory

    if ( NULL != recv_buf  ) std::free( recv_buf );
    if ( NULL != recv_len  ) std::free( recv_len );
    if ( NULL != recv_disp ) std::free( recv_disp );
  }

  comm_time += MPI_Wtime() - start_time ;

  return -1 == result ? (mpi_filebuf *) NULL : this ;
}

/*--------------------------------------------------------------------*/

int mpi_filebuf::sync()
{
  // The root processor will push to file, all others ignore

  if ( NULL != comm_root_fp ) {

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    int    cur_len = pptr() - cur_buf ;

    if ( 0 < cur_len ) std::fwrite(cur_buf,1,cur_len,comm_root_fp);

    std::fflush( comm_root_fp );

    setp( comm_buffer , epptr() );
  }

  return 0 ;
}


std::streambuf * mpi_filebuf::setbuf( char * s , std::streamsize n )
{
  return this ;
}
