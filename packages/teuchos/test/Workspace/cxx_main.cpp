// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include <valarray>

#include "Teuchos_Workspace.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Version.hpp"

/** \brief This class implements a simple (useless) transformation that requires
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

  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    // Read options from the commandline
    CommandLineProcessor  clp(false); // Don't throw exceptions

    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );

    double rel_proc_speed = 1e-5; // Should 
    clp.setOption( "rel-proc-speed", &rel_proc_speed, "Relative processor speed (try around 1.0 for timing)." );

    int size = 1;
    clp.setOption( "size", &size, "Size of memory blocks created." );

    bool allocate_workspace = true;
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

    if (verbose)
      std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

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
      Teuchos::print_memory_usage_stats(Teuchos::get_default_workspace_store().get(),std::cout);
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
			std::cerr << "*** Caught standard std::exception : " << excpt.what() << std::endl;
		return 1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown std::exception\n";
		return 1;
	}
  
	return 0;
  
}
