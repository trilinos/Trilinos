// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <valarray>

#include "Teuchos_Workspace.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"

///
/** This class implements a simple (useless) transformation that requires
 * workspace.
 *
 * This class creates workspace using one of four approaches:
 * <ul>
 * <li> Using raw calls to new [] and delete []
 * <li> Using std::vector
 * <li> Using std::valarray
 * <li> Using Teuchos::Workspace
 * </ul>
 */
class Transformer {
  Teuchos::WorkspaceStore *wss_;
  void transform( const int size, double a[], double b[] ) {
    b[0] = a[0];
    for( int k = 1; k < size; ++k )  b[k]  = a[k]+a[k-1];
    for( int k = 0; k < size; ++k )  a[k]  = a[k]-b[k];
  }
public:
  Transformer() : wss_(Teuchos::get_default_workspace_store().get()) {}
  void transformRaw( const int size, double a[] ) {
    double *b = new double[size]; // Should not call constructors!
    transform( size, a, b );
    delete [] b;
  }
  void transformVector( const int size, double a[] ) {
    std::vector<double> b(size); // Should call constructors!
    transform( size, a, &b[0] );
  }
  void transformValarray( const int size, double a[] ) {
    std::valarray<double> b(size); // Should not call constructors!
    transform( size, a, &b[0] );
  }
  void transformWorkspace( const int size, double a[] ) {
    Teuchos::Workspace<double> b(wss_,size,false); // Does not call constructors!
    transform( size, a, &b[0] );
  }
};

int main( int argc, char* argv[] )
{
  using Teuchos::CommandLineProcessor;

  bool verbose = false;

  try {

    // Read options from the commandline
    double rel_proc_speed   = 1.0;
    int  size               = 1;
    bool allocate_workspace = true;
    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "rel-proc-speed", &rel_proc_speed, "Relative processor speed." );
    clp.setOption( "size", &size, "Size of memory blocks created." );
    clp.setOption( "allocate-workspace", "no-allocate-workspace", &allocate_workspace, "Preallocate workspace or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    // Determine how many loops to do to get good timings
    const long int
      default_num_loops = int( 100000000 * rel_proc_speed ),
      num_loops         = int( default_num_loops / ( size + 100 ) );

    // Allocate workspace
    if( allocate_workspace )
      Teuchos::set_default_workspace_store(
        Teuchos::rcp(new Teuchos::WorkspaceStoreInitializeable(10*size))
        );

		Teuchos::Time timer("");

    if(verbose) std::cout
      << "\n************************************************************************************"
      << "\n*** Testing and timing Teuchos::Workspace and other methods for temporary memory ***"
      << "\n************************************************************************************\n";

    if(verbose) std::cout
      << "\nMemory block size    = " << size
      << "\nNumber of call loops = " << num_loops
      << std::endl;

    Transformer t;
    std::vector<double> a(size);

    if(verbose) std::cout << "\nTiming raw new and delete for temporaries ...\n";
    std::fill_n( &a[0], size, 1.0 );
		timer.start(true);
    for( int k = 0; k < num_loops; ++k ) t.transformRaw(size,&a[0]);
    timer.stop();
    const double raw_time = timer.totalElapsedTime();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";

    if(verbose) std::cout << "\nTiming std::vector for temporaries ...\n";
    std::fill_n( &a[0], size, 1.0 );
		timer.start(true);
    for( int k = 0; k < num_loops; ++k ) t.transformVector(size,&a[0]);
    timer.stop();
    const double vector_time = timer.totalElapsedTime();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";

    if(verbose) std::cout << "\nTiming std::valarray for temporaries ...\n";
    std::fill_n( &a[0], size, 1.0 );
		timer.start(true);
    for( int k = 0; k < num_loops; ++k ) t.transformValarray(size,&a[0]);
    timer.stop();
    const double valarray_time = timer.totalElapsedTime();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";

    if(verbose) std::cout << "\nTiming Teuchos::Workspace for temporaries ...\n";
    std::fill_n( &a[0], size, 1.0 );
		timer.start(true);
    for( int k = 0; k < num_loops; ++k ) t.transformWorkspace(size,&a[0]);
    timer.stop();
    const double workspace_time = timer.totalElapsedTime();
		if(verbose) std::cout << "  time = " << timer.totalElapsedTime() << " sec\n";

    if(verbose) {
      std::cout
        << "\nRelative time (lower is better):"
        << "\n   raw new/delete      = " << (raw_time/workspace_time)
        << "\n   std::vector         = " << (vector_time/workspace_time)
        << "\n   std::valarray       = " << (valarray_time/workspace_time)
        << "\n   Teuchos::Workspace  = " << (workspace_time/workspace_time)
        << std::endl << std::endl;
    }

  }
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
		return 1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception\n";
		return 1;
	}
  
	return 0;
  
}
