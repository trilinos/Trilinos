// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef RTOp_USE_MPI
#  include "Epetra_MpiComm.h"
#endif

//
// Some helper functions
//

namespace {

void print_performance_stats(
  const int        num_time_samples
  ,const double    raw_epetra_time
  ,const double    thyra_wrapped_time
  ,bool            verbose
  ,std::ostream    &out
  )
{
  if(verbose)
    out
      << "\nAverage times (out of " << num_time_samples << " samples):\n"
      << "  Raw Epetra              = " << (raw_epetra_time/num_time_samples) << std::endl
      << "  Thyra Wrapped Epetra    = " << (thyra_wrapped_time/num_time_samples) << std::endl
      << "\nRelative performance of Thyra wrapped verses raw Epetra:\n"
      << "  ( raw epetra time / thyra wrapped time ) = ( " << raw_epetra_time << " / " << thyra_wrapped_time << " ) = "
      << (raw_epetra_time/thyra_wrapped_time) << std::endl;
}

inline
double sum( const Epetra_MultiVector &ev )
{
  std::vector<double> meanValues(ev.NumVectors());
  ev.MeanValue(&meanValues[0]);
  double sum = 0;
  for( int i = 0; i < ev.NumVectors(); ++i ) sum += meanValues[i];
  return (ev.Map().NumGlobalElements()*sum);
}

} // namespace

/* Testing program for Thyra/Epetra adpaters.
 *
 * This testing program shows how you can easily mix and match
 * different implementations of vectors and multi-vectors for serial
 * and SPMD MPI implementations.  This code is worth study to show how
 * this is done.
 *
 * Note that the tests performed do not prove that the Epetra adapters
 * (or Epetra itself) perform correctly as only a few post conditions
 * are checked.  Because of the simple nature of these computations it
 * would be possible to put in more exactly component-wise tests if
 * that is needed in the future.
 */
int main( int argc, char* argv[] )
{

  using std::endl;

  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  using Teuchos::dyn_cast;
  using Teuchos::CommandLineProcessor;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::arcpFromArray;
  using Teuchos::rcp;
  using Teuchos::rcp_static_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::testRelErr;

  using Thyra::passfail;
  using Thyra::NOTRANS;
  using Thyra::TRANS;
  using Thyra::apply;
  using Thyra::create_VectorSpace;
  using Thyra::create_Vector;
  using Thyra::create_MultiVector;
  
  bool verbose = true;
  bool dumpAll = false;
  bool success = true;
  bool result;

  int procRank = 0;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {

    Teuchos::Time timer("");

    //
    // Read options from the commandline
    //

    int     local_dim            = 1000;
    int     num_mv_cols          = 4;
    double  max_rel_err          = 1e-13;
    double  max_rel_warn         = 1e-15;
    double  scalar               = 1.5;
    double  max_flop_rate        = 1.0; // Turn off by default!
#ifdef RTOp_USE_MPI
    bool    useMPI               = true;
#endif
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
    clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );
    clp.setOption( "num-mv-cols", &num_mv_cols, "Number columns in each multi-vector (>=4)." );
    clp.setOption( "max-rel-err-tol", &max_rel_err, "Maximum relative error tolerance for tests." );
    clp.setOption( "max-rel-warn-tol", &max_rel_warn, "Maximum relative warning tolerance for tests." );
    clp.setOption( "scalar", &scalar, "A scalar used in all computations." );
    clp.setOption( "max-flop-rate", &max_flop_rate, "Approx flop rate used for loop timing." );
#ifdef RTOp_USE_MPI
    clp.setOption( "use-mpi", "no-use-mpi", &useMPI, "Actually use MPI or just run independent serial programs." );
#endif
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPTION(
      num_mv_cols < 4, std::logic_error
      ,"Error, num-mv-cols must be >= 4!"
      );

    //
    // Get basic MPI info
    //

#ifdef RTOp_USE_MPI
    MPI_Comm mpiComm;
    int numProc;
    if(useMPI) {
      mpiComm = MPI_COMM_WORLD;
      MPI_Comm_size( mpiComm, &numProc );
      MPI_Comm_rank( mpiComm, &procRank );
    }
    else {
      mpiComm = MPI_COMM_NULL;
      numProc = 1;
      procRank = 0;
    }
#endif

    if(verbose)
      *out
        << "\n***"
        << "\n*** (A) Creating two vector spaces (an Epetra-based and a non-Epetra-based)"
        << "\n***\n";

    //
    // Create two different vector spaces (one Epetra and one
    // non-Epetra) that should be compatible.
    //
    RCP<const Epetra_Comm> epetra_comm;
    RCP<const Epetra_Map> epetra_map;
    RCP<const Thyra::VectorSpaceBase<Scalar> > epetra_vs;
    RCP<const Thyra::VectorSpaceBase<Scalar> > non_epetra_vs;
#ifdef RTOp_USE_MPI
    if(useMPI) {
      //
      // Create parallel vector spaces with compatible maps
      //
      // Epetra vector space
      if(verbose) *out << "\nCreating vector space using Epetra_MpiComm ...\n";
      epetra_comm = rcp(new Epetra_MpiComm(mpiComm));
      epetra_map = rcp(new Epetra_Map(-1,local_dim,0,*epetra_comm));
      epetra_vs = Thyra::create_VectorSpace(epetra_map);
      // Non-Epetra vector space
      if(verbose) *out << "\nCreating Thyra::DefaultSpmdVectorSpace ...\n";
      non_epetra_vs = rcp(
        new Thyra::DefaultSpmdVectorSpace<Scalar>(
          Thyra::create_Comm(epetra_comm)
          ,local_dim,-1
          )
        );
    }
    else {
#endif
      //
      // Create serial vector spaces (i.e. VectorSpaceBase::isInCore()==true)
      //
      // Epetra vector space
      if(verbose) *out << "\nCreating vector space using Epetra_SerialComm ...\n";
      epetra_comm = rcp(new Epetra_SerialComm);
      epetra_map = rcp(new Epetra_LocalMap(local_dim,0,*epetra_comm));
      epetra_vs = Thyra::create_VectorSpace(epetra_map);
      // Non-Epetra vector space
      if(verbose) *out << "\nCreating Thyra::DefaultSpmdVectorSpace ...\n";
      non_epetra_vs = Thyra::defaultSpmdVectorSpace<Scalar>(local_dim);
#ifdef RTOp_USE_MPI
    }
#endif // end create vector spacdes [Doxygen looks for this!]

#ifdef RTOp_USE_MPI
    const int global_dim = local_dim * numProc;
#else
    const int global_dim = local_dim;
#endif

    if(verbose)
      *out
        << "\nscalar              = " << scalar
        << "\nlocal_dim           = " << local_dim
        << "\nglobal_dim          = " << global_dim
        << "\nnum_mv_cols         = " << num_mv_cols
        << "\nepetra_vs.dim()     = " << epetra_vs->dim()
        << "\nnon_epetra_vs.dim() = " << non_epetra_vs->dim()
        << std::endl;

    //
    // Create vectors and multi-vectors from each vector space
    //

    RCP<Thyra::VectorBase<Scalar> >
      ev1 = createMember(epetra_vs),
      ev2 = createMember(epetra_vs);
    RCP<Thyra::VectorBase<Scalar> >
      nev1 = createMember(non_epetra_vs),
      nev2 = createMember(non_epetra_vs);

    RCP<Thyra::MultiVectorBase<Scalar> >
      eV1 = createMembers(epetra_vs,num_mv_cols),
      eV2 = createMembers(epetra_vs,num_mv_cols);
    RCP<Thyra::MultiVectorBase<Scalar> >
      neV1 = createMembers(non_epetra_vs,num_mv_cols),
      neV2 = createMembers(non_epetra_vs,num_mv_cols);

    if(verbose)
      *out
        << "\n***"
        << "\n*** (B) Testing Epetra and non-Epetra Thyra wrapped objects"
        << "\n***\n";

    //
    // Check for compatibility of the vector and Multi-vectors
    // w.r.t. RTOps
    //

    if(verbose) *out << "\n*** (B.1) Testing individual vector/multi-vector RTOps\n";

    Thyra::assign( ev1.ptr(), 0.0 );
    Thyra::assign( ev2.ptr(), scalar );
    Thyra::assign( nev1.ptr(), 0.0 );
    Thyra::assign( nev2.ptr(), scalar );
    Thyra::assign( eV1.ptr(), 0.0 );
    Thyra::assign( eV2.ptr(), scalar );
    Thyra::assign( neV1.ptr(), 0.0 );
    Thyra::assign( neV2.ptr(), scalar );

    Scalar
      ev1_nrm = Thyra::norm_1(*ev1),
      ev2_nrm = Thyra::norm_1(*ev2),
      eV1_nrm = Thyra::norm_1(*eV1),
      eV2_nrm = Thyra::norm_1(*eV2),
      nev1_nrm = Thyra::norm_1(*nev1),
      nev2_nrm = Thyra::norm_1(*nev2),
      neV1_nrm = Thyra::norm_1(*neV1),
      neV2_nrm = Thyra::norm_1(*neV2);

    const std::string s1_n = "fabs(scalar)*global_dim";
    const Scalar s1 = fabs(scalar)*global_dim;
    
    if(!testRelErr("Thyra::norm_1(ev1)",ev1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nev1 =\n" << *ev1;
    if(!testRelErr("Thyra::norm_1(ev2)",ev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nev2 =\n" << *ev2;
    if(!testRelErr("Thyra::norm_1(nev1)",nev1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nnev2 =\n" << *ev1;
    if(!testRelErr("Thyra::norm_1(nev2)",nev2_nrm,s1_n,s1,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nnev2 =\n" << *nev2;
    if(!testRelErr("Thyra::norm_1(eV1)",eV1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV1 =\n" << *eV1;
    if(!testRelErr("Thyra::norm_1(eV2)",eV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV2 =\n" << *eV2;
    if(!testRelErr("Thyra::norm_1(neV1)",neV1_nrm,"0",Scalar(0),"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV1 =\n" << *neV1;
    if(!testRelErr("Thyra::norm_1(neV2)",neV2_nrm,s1_n,s1,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV2 =\n" << *neV2;

    if(verbose) *out << "\n*** (B.2) Test RTOps with two or more arguments\n";

    if(verbose) *out << "\nPerforming ev1 = ev2 ...\n";
    timer.start(true);
     Thyra::assign( ev1.ptr(), *ev2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ev1)",Thyra::norm_1(*ev1),"Thyra::norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nev1 =\n" << *ev1;

    if(verbose) *out << "\nPerforming eV1 = eV2 ...\n";
    timer.start(true);
     Thyra::assign( eV1.ptr(), *eV2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eV1)",Thyra::norm_1(*eV1),"Thyra::norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV1 =\n" << *eV1;

    if(verbose) *out << "\nPerforming ev1 = nev2 ...\n";
    timer.start(true);
     Thyra::assign( ev1.ptr(), *nev2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ev1)",Thyra::norm_1(*ev1),"Thyra::norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nev1 =\n" << *ev1;

    if(verbose) *out << "\nPerforming nev1 = ev2 ...\n";
    timer.start(true);
     Thyra::assign( nev1.ptr(), *ev2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(nev1)",Thyra::norm_1(*nev1),"Thyra::norm_1(ev2)",ev2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nnev1 =\n" << *nev1;

    if(verbose) *out << "\nPerforming nev1 = nev2 ...\n";
    timer.start(true);
     Thyra::assign( nev1.ptr(), *nev2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(nev1)",Thyra::norm_1(*nev1),"Thyra::norm_1(nev2)",nev2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nnev1 =\n" << *nev1;

    if(verbose) *out << "\nPerforming eV1 = neV2 ...\n";
    timer.start(true);
     Thyra::assign( eV1.ptr(), *neV2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eV1)",Thyra::norm_1(*eV1),"Thyra::norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV1 =\n" << *eV1;

    if(verbose) *out << "\nPerforming neV1 = eV2 ...\n";
    timer.start(true);
     Thyra::assign( neV1.ptr(), *eV2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neV1)",Thyra::norm_1(*neV1),"Thyra::norm_1(eV2)",eV2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV1 =\n" << *neV1;

    if(verbose) *out << "\nPerforming neV1 = neV2 ...\n";
    timer.start(true);
     Thyra::assign( neV1.ptr(), *neV2 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neV1)",Thyra::norm_1(*neV1),"Thyra::norm_1(neV2)",neV2_nrm,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV1 =\n" << *neV1;

    Thyra::LinearOpTester<Scalar> linearOpTester;
    linearOpTester.set_all_warning_tol(max_rel_warn);
    linearOpTester.set_all_error_tol(max_rel_err);
    linearOpTester.dump_all(dumpAll);

    Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
    linearOpWithSolveTester.set_all_solve_tol(max_rel_err);
    linearOpWithSolveTester.set_all_slack_error_tol(max_rel_err);
    linearOpWithSolveTester.set_all_slack_warning_tol(max_rel_warn);
    linearOpWithSolveTester.dump_all(dumpAll);


    if(verbose) *out << "\n*** (B.3) Test Vector linear operator interface\n";

    if(verbose) *out << "\nChecking *out linear operator interface of ev1 ...\n";
    result = linearOpTester.check(*ev1,out.ptr());
    if(!result) success = false;

    if(verbose) *out << "\nChecking *out linear operator interface of nev1 ...\n";
    result = linearOpTester.check(*nev1,out.ptr());
    if(!result) success = false;


    if(verbose) *out << "\n*** (B.4) Test MultiVector linear operator interface\n";

    if(verbose) *out << "\nChecking *out linear operator interface of eV1 ...\n";
    result = linearOpTester.check(*eV1,out.ptr());
    if(!result) success = false;

    if(verbose) *out << "\nChecking *out linear operator interface of neV1 ...\n";
    result = linearOpTester.check(*neV1,out.ptr());
    if(!result) success = false;

    const std::string s2_n = "scalar^2*global_dim*num_mv_cols";
    const Scalar s2 = scalar*scalar*global_dim*num_mv_cols;

    RCP<Thyra::MultiVectorBase<Scalar> >
      T = createMembers(eV1->domain(),num_mv_cols);


    if(verbose) *out << "\n*** (B.5) Test MultiVector::apply(...)\n";

    if(verbose) *out << "\nPerforming eV1'*eV2 ...\n";
    timer.start(true);
    apply( *eV1, TRANS, *eV2, T.ptr() );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eV1'*eV2)",Thyra::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV1'*eV2 =\n" << *T;

    if(verbose) *out << "\nPerforming neV1'*eV2 ...\n";
    timer.start(true);
    apply( *neV1, TRANS, *eV2, T.ptr() );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neV1'*eV2)",Thyra::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV1'*eV2 =\n" << *T;

    if(verbose) *out << "\nPerforming eV1'*neV2 ...\n";
    timer.start(true);
    apply( *eV1, TRANS, *neV2, T.ptr() );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eV1'*neV2)",Thyra::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\neV1'*neV2 =\n" << *T;

    if(verbose) *out << "\nPerforming neV1'*neV2 ...\n";
    timer.start(true);
    apply( *neV1, TRANS, *neV2, T.ptr() );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neV1'*neV2)",Thyra::norm_1(*T),s2_n,s2,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
    if(verbose && dumpAll) *out << "\nneV1'*neV2 =\n" << *T;


    if(verbose) *out << "\n*** (B.6) Creating a diagonal Epetra_Operator Op\n";

    RCP<Epetra_Operator>  epetra_op;

    {
      // Create a diagonal matrix with scalar on the diagonal
      RCP<Epetra_CrsMatrix>
        epetra_mat = rcp(new Epetra_CrsMatrix(::Copy,*epetra_map,1));
      Scalar values[1] = { scalar };
      int indices[1];
      const int IB = epetra_map->IndexBase(), offset = procRank*local_dim;
      for( int k = 0; k < local_dim; ++k ) {
        indices[0] = offset + k + IB;  // global column
        epetra_mat->InsertGlobalValues(
          offset + k + IB     // GlobalRow
          ,1                  // NumEntries
          ,values             // Values
          ,indices            // Indices
          );
      }
      epetra_mat->FillComplete();
      epetra_op = epetra_mat;
    } // end epetra_op

    RCP<const Thyra::LinearOpBase<Scalar> >
      Op = Thyra::epetraLinearOp(epetra_op);

    if(verbose && dumpAll) *out << "\nOp=\n" << *Op;


    if(verbose) *out << "\n*** (B.6b) Going through partial then full initialization of EpetraLinearOp ...\n";

    {
      if(verbose) *out
        << "\nChecking isFullyUninitialized(*nonconstEpetraLinearOp())==true : ";
      RCP<Thyra::EpetraLinearOp>
        thyraEpetraOp = Thyra::nonconstEpetraLinearOp();
      result = isFullyUninitialized(*thyraEpetraOp);
      if(!result) success = false;
      if(verbose) *out << Thyra::passfail(result) << endl;
    }

    {

      if(verbose) *out
        << "\nthyraEpetraOp = partialNonconstEpetraLinearOp(...)\n";
      RCP<Thyra::EpetraLinearOp> thyraEpetraOp =
        Thyra::partialNonconstEpetraLinearOp(
          epetra_vs, epetra_vs, epetra_op
          );

      if(verbose) *out
        << "\nChecking isPartiallyInitialized(*thyraEpetraOp)==true : ";
      result = isPartiallyInitialized(*thyraEpetraOp);
      if(!result) success = false;
      if(verbose) *out << Thyra::passfail(result) << endl;

      if(verbose) *out
        << "\nthyraEpetraOp->setFullyInitialized(true)\n";
      thyraEpetraOp->setFullyInitialized(true);

      if(verbose) *out
        << "\nChecking isFullyInitialized(*thyraEpetraOp)==true : ";
      result = isFullyInitialized(*thyraEpetraOp);
      if(!result) success = false;
      if(verbose) *out << Thyra::passfail(result) << endl;

    }


    if(verbose) *out << "\n*** (B.7) Test EpetraLinearOp linear operator interface\n";

    if(verbose) *out << "\nChecking *out linear operator interface of Op ...\n";
    result = linearOpTester.check(*Op,out.ptr());
    if(!result) success = false;

    RCP<Thyra::VectorBase<Scalar> >
      ey  = createMember(epetra_vs);
    RCP<Thyra::MultiVectorBase<Scalar> >
      eY  = createMembers(epetra_vs,num_mv_cols);
    RCP<Thyra::VectorBase<Scalar> >
      ney = createMember(non_epetra_vs);
    RCP<Thyra::MultiVectorBase<Scalar> >
      neY = createMembers(non_epetra_vs,num_mv_cols);


    if(verbose) *out << "\n*** (B.8) Mix and match vector and Multi-vectors with Epetra opeator\n";

    const std::string s3_n = "2*scalar^2*global_dim";
    const Scalar s3 = 2*scalar*scalar*global_dim;
    
    if(verbose) *out << "\nPerforming ey = 2*Op*ev1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *ev1, ey.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ey)",Thyra::norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming eY = 2*Op*eV1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1, eY.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eY)",Thyra::norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming ney = 2*Op*ev1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *ev1, ney.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ney)",Thyra::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming neY = 2*Op*eV1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1, neY.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neY)",Thyra::norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming ey = 2*Op*nev1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *nev1, ey.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ey)",Thyra::norm_1(*ey),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming eY = 2*Op*neV1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *neV1, eY.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eY)",Thyra::norm_1(*eY),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming ney = 2*Op*nev1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *nev1, ney.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ney)",Thyra::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming ney = 2*Op*nev1 through MultiVector interface ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, static_cast<const Thyra::MultiVectorBase<Scalar>&>(*nev1), Ptr<Thyra::MultiVectorBase<Scalar> >(ney.ptr()), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(ney)",Thyra::norm_1(*ney),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming neY = 2*Op*neV1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *neV1, neY.ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(neY)",Thyra::norm_1(*neY),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;


    if(verbose) *out << "\n*** (B.9) Testing Multi-vector views with Epetra operator\n";

    const Thyra::Range1D col_rng(0,1);
    const Teuchos::Tuple<int, 2> cols = Teuchos::tuple<int>(2, 3);

    RCP<const Thyra::MultiVectorBase<Scalar> >
      eV1_v1  = rcp_static_cast<const Thyra::MultiVectorBase<Scalar> >(eV1)->subView(col_rng),
      eV1_v2  = rcp_static_cast<const Thyra::MultiVectorBase<Scalar> >(eV1)->subView(cols);
    RCP<const Thyra::MultiVectorBase<Scalar> >
      neV1_v1  = rcp_static_cast<const Thyra::MultiVectorBase<Scalar> >(neV1)->subView(col_rng),
      neV1_v2  = rcp_static_cast<const Thyra::MultiVectorBase<Scalar> >(neV1)->subView(cols);
    if(verbose && dumpAll) {
      *out << "\neV1_v1=\n" << *eV1_v1;
      *out << "\neV1_v2=\n" << *eV1_v2;
      *out << "\nneV1_v1=\n" << *neV1_v1;
      *out << "\nneV1_v2=\n" << *neV1_v2;
    }

    if(verbose) *out << "\nPerforming eY_v1 = 2*Op*eV1_v1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1_v1, eY->subView(col_rng).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(verbose && dumpAll) *out << "\neV_v1=\n" << *eY->subView(col_rng);
    if(!testRelErr("Thyra::norm_1(eY_v1)",Thyra::norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming eY_v2 = 2*Op*eV1_v2 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1_v2, eY->subView(cols).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(verbose && dumpAll) *out << "\neV_v2=\n" << *eY->subView(cols);
    if(!testRelErr("Thyra::norm_1(eY_v2)",Thyra::norm_1(*eY->subView(cols)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming neY_v1 = 2*Op*eV1_v1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1_v1, neY->subView(col_rng).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(verbose && dumpAll) *out << "\neV_v1=\n" << *eY->subView(col_rng);
    if(!testRelErr("Thyra::norm_1(neY_v1)",Thyra::norm_1(*neY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming eY_v1 = 2*Op*neV1_v1 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *neV1_v1, eY->subView(col_rng).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(!testRelErr("Thyra::norm_1(eY_v1)",Thyra::norm_1(*eY->subView(col_rng)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming neY_v2 = 2*Op*eV1_v2 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *eV1_v2, neY->subView(cols).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(verbose && dumpAll) *out << "\neV_v2=\n" << *eY->subView(cols);
    if(!testRelErr("Thyra::norm_1(neY_v2)",Thyra::norm_1(*neY->subView(cols)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;

    if(verbose) *out << "\nPerforming eY_v2 = 2*Op*neV1_v2 ...\n";
    timer.start(true);
    apply( *Op, NOTRANS, *neV1_v2, eY->subView(cols).ptr(), 2.0 );
    timer.stop();
    if(verbose) *out << "  time = " << timer.totalElapsedTime() << " sec\n";
    if(verbose && dumpAll) *out << "\neV_v2=\n" << *eY->subView(cols);
    if(!testRelErr("Thyra::norm_1(eY_v2)",Thyra::norm_1(*eY->subView(cols)),s3_n,s3,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;


    if(verbose) *out << "\n*** (B.10) Testing Vector and MultiVector view creation functions\n";

    {

      const std::string s_n = "fabs(scalar)*num_mv_cols";
      const Scalar s = fabs(scalar)*num_mv_cols;

      Array<Scalar> t_raw_values( num_mv_cols );
      RTOpPack::SubVectorView<Scalar> t_raw( 0, num_mv_cols,
        arcpFromArray(t_raw_values), 1 );

      std::fill_n( t_raw_values.begin(), t_raw_values.size(), ST::zero() );
      Thyra::assign( createMemberView(T->range(),t_raw).ptr(), scalar );
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > t_view = createMemberView(T->range(),static_cast<RTOpPack::ConstSubVectorView<Scalar>&>(t_raw));
      Scalar t_nrm = Thyra::norm_1(*t_view);
      if(!testRelErr("Thyra::norm_1(t_view)",t_nrm,s_n,s,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
      if(verbose && dumpAll) *out << "\nt_view =\n" << *t_view;

/*
#ifndef SUN_CXX // The sun compiler Forte Developer 5.4 does not destory temporaries properly and this does not work
      std::fill_n( t_raw_values.begin(), t_raw_values.size(), ST::zero() );
      Thyra::assign( T->range().ptr()->Thyra::VectorSpaceBase<Scalar>::createMemberView(t_raw), scalar );
      t_view = T->range()->Thyra::VectorSpaceBase<Scalar>::createMemberView(static_cast<RTOpPack::ConstSubVectorView<Scalar>&>(t_raw));
      t_nrm = Thyra::norm_1(*t_view);
      if(!testRelErr("Thyra::norm_1(t_view)",t_nrm,s_n,s,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
      if(verbose && dumpAll) *out << "\nt_view =\n" << *t_view;
#endif
*/

      Array<Scalar> T_raw_values( num_mv_cols * num_mv_cols );
      RTOpPack::SubMultiVectorView<Scalar> T_raw( 0, num_mv_cols, 0, num_mv_cols,
        arcpFromArray(T_raw_values), num_mv_cols );

      std::fill_n( T_raw_values.begin(), T_raw_values.size(), ST::zero() );
      Thyra::assign( createMembersView(T->range(),T_raw).ptr(), scalar );
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
        T_view = createMembersView(T->range(),static_cast<RTOpPack::ConstSubMultiVectorView<Scalar>&>(T_raw));
      Scalar T_nrm = Thyra::norm_1(*T_view);
      if(!testRelErr("Thyra::norm_1(T_view)",T_nrm,s_n,s,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
      if(verbose && dumpAll) *out << "\nT_view =\n" << *T_view;

/*
#ifndef SUN_CXX // The sun compiler Forte Developer 5.4 does not destory temporaries properly and this does not work
      std::fill_n( T_raw_values.begin(), T_raw_values.size(), ST::zero() );
      Thyra::assign( T->range().ptr()->Thyra::VectorSpaceBase<Scalar>::createMembersView(T_raw), scalar );
      T_view = T->range()->Thyra::VectorSpaceBase<Scalar>::createMembersView(static_cast<RTOpPack::ConstSubMultiVectorView<Scalar>&>(T_raw));
      T_nrm = Thyra::norm_1(*T_view);
      if(!testRelErr("Thyra::norm_1(T_view)",T_nrm,s_n,s,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())) success=false;
      if(verbose && dumpAll) *out << "\nT_view =\n" << *T_view;
#endif
*/

    }


    if(verbose) *out << "\n*** (B.11) Testing Epetra_Vector and Epetra_MultiVector wrappers\n";

    {

      Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<Scalar> >
        mpi_vs = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar> >(epetra_vs,true);

      if(verbose) *out << "\nCreating temporary Epetra_Vector et1 and Epetra_MultiVector eT1 objects ...\n";
      Teuchos::RCP<Epetra_Vector>
        et1 = Teuchos::rcp(new Epetra_Vector(*epetra_map));
      Teuchos::RCP<Epetra_MultiVector>
        eT1 = Teuchos::rcp(new Epetra_MultiVector(*epetra_map,num_mv_cols));

      if(verbose) *out << "\nCreating non-const VectorBase t1 and MultiVectorBase T1 objects from et1 and eT2 ...\n";
      Teuchos::RCP<Thyra::VectorBase<Scalar> >
        t1 = create_Vector(et1,mpi_vs);
      Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
        T1 = create_MultiVector(eT1,mpi_vs);

      if(verbose) *out << "\nPerforming t1 = ev1 ...\n";
      assign( t1.ptr(), *ev1 );
      if(!testRelErr(
           "sum(t1)",Thyra::sum(*t1),"sum(ev1)",sum(*ev1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(verbose) *out << "\nPerforming T1 = eV1 ...\n";
      assign( T1.ptr(), *eV1 );
      if(!testRelErr(
           "norm_1(T1)",Thyra::norm_1(*T1),"norm_1(eV1)",norm_1(*eV1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(verbose) *out << "\nChecking that t1 and T1 yield the same objects as et1 and eT2 ...\n";
  
      Teuchos::RCP<Epetra_Vector>
        et1_v = get_Epetra_Vector(*epetra_map,t1);
      result = et1_v.get() == et1.get();
      if(verbose) *out << "\net1_v.get() = " << et1_v.get() << " == et1.get() = " << et1.get() << " : " << passfail(result) << endl;
      if(!result) success = false;

      Teuchos::RCP<Epetra_MultiVector>
        eT1_v = get_Epetra_MultiVector(*epetra_map,T1);
      result = eT1_v.get() == eT1.get();
      if(verbose) *out << "\neT1_v.get() = " << eT1_v.get() << " == eT1.get() = " << eT1.get() << " : " << passfail(result) << endl;
      if(!result) success = false;

      if(verbose) *out << "\nCreating const VectorBase ct1 and MultiVectorBase cT1 objects from et1 and eT2 ...\n";
      Teuchos::RCP<const Thyra::VectorBase<Scalar> >
        ct1 = create_Vector(Teuchos::rcp_implicit_cast<const Epetra_Vector>(et1),mpi_vs);
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
        cT1 = create_MultiVector(Teuchos::rcp_implicit_cast<const Epetra_MultiVector>(eT1),mpi_vs);

      if(!testRelErr(
           "sum(ct1)",Thyra::sum(*ct1),"sum(ev1)",sum(*ev1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(!testRelErr(
           "norm_1(cT1)",Thyra::norm_1(*cT1),"norm_1(eV1)",norm_1(*eV1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(verbose) *out << "\nChecking that ct1 and cT1 yield the same objects as et1 and eT2 ...\n";
  
      Teuchos::RCP<const Epetra_Vector>
        cet1_v = get_Epetra_Vector(*epetra_map,ct1);
      result = cet1_v.get() == et1.get();
      if(verbose) *out << "\ncet1_v.get() = " << cet1_v.get() << " == et1.get() = " << et1.get() << " : " << passfail(result) << endl;
      if(!result) success = false;

      Teuchos::RCP<const Epetra_MultiVector>
        ceT1_v = get_Epetra_MultiVector(*epetra_map,cT1);
      result = ceT1_v.get() == eT1.get();
      if(verbose) *out << "\nceT1_v.get() = " << ceT1_v.get() << " == eT1.get() = " << eT1.get() << " : " << passfail(result) << endl;
      if(!result) success = false;

      if(verbose) *out << "\nCreating non-const Epetra_Vector ett1 and Epetra_MultiVector etT1 objects from clones of t1 and T2 ...\n";
      Teuchos::RCP<Epetra_Vector>
        ett1 = get_Epetra_Vector(*epetra_map,t1->clone_v());
      Teuchos::RCP<Epetra_MultiVector>
        etT1 = get_Epetra_MultiVector(*epetra_map,T1->clone_mv());

      if(verbose) *out << "\nChecking that ett1 and etT1 yield objects with the same value as et1 and eT2 ...\n";

      if(!testRelErr(
           "sum(ett1)",sum(*ett1),"sum(et1)",sum(*et1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(!testRelErr(
           "sum(etT1)",sum(*etT1),"sum(eT1)",sum(*eT1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(verbose) *out << "\nCreating const Epetra_Vector cett1 and Epetra_MultiVector cetT1 objects from clones of t1 and T2 ...\n";
      Teuchos::RCP<const Epetra_Vector>
        cett1 = get_Epetra_Vector(*epetra_map,Teuchos::rcp_implicit_cast<const Thyra::VectorBase<Scalar> >(t1->clone_v()));
      Teuchos::RCP<const Epetra_MultiVector>
        cetT1 = get_Epetra_MultiVector(*epetra_map,Teuchos::rcp_implicit_cast<const Thyra::MultiVectorBase<Scalar> >(T1->clone_mv()));

      if(verbose) *out << "\nChecking that cett1 and cetT1 yield objects with the same value as et1 and eT2 ...\n";

      if(!testRelErr(
           "sum(cett1)",sum(*cett1),"sum(et1)",sum(*et1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

      if(!testRelErr(
           "sum(cetT1)",sum(*cetT1),"sum(eT1)",sum(*eT1)
           ,"max_rel_err",max_rel_err,"max_rel_warn",max_rel_warn,out.ptr())
        ) success=false;

    }


    if(verbose) *out << "\n*** (B.12) Test DiagonalEpetraLinearOpWithSolveFactory \n";

    {

      if(verbose) *out << "\nUsing DiagonalEpetraLinearOpWithSolveFactory to create diagLOWS from Op ...\n";
      
      Thyra::DiagonalEpetraLinearOpWithSolveFactory diagLOWSFactory;

      Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
        diagLOWS = Thyra::linearOpWithSolve<double>(diagLOWSFactory,Op);

      if(verbose) *out << "\nTesting LinearOpBase interface of diagLOWS ...\n";

      result = linearOpTester.check(*diagLOWS, out.ptr());
      if(!result) success = false;

      if(verbose) *out << "\nTesting LinearOpWithSolveBase interface of diagLOWS ...\n";

      result = linearOpWithSolveTester.check(*diagLOWS, &*out);
      if(!result) success = false;
    
    }


    if(verbose)
      *out
        << "\n***"
        << "\n*** (C) Comparing the speed of Thyra adapted Epetra objects verses raw Epetra objects"
        << "\n***\n";

    //
    // Setup the number of timing loops to get good timings
    //
    // Here we try to shoot for timing ab*out 1 second's worth of
    // computations and adjust the number of evaluation loops
    // accordingly.  Let X be the approximate number of flops per
    // loop (or per evaluation).  We then compute the number of
    // loops as:
    //
    // 1.0 sec |  num CPU flops |   1 loop  |
    // --------|----------------|-----------|
    //         |       sec      |   X flops |
    //
    // This just comes *out to be:
    //
    //   num_time_loops_X =  max_flop_rate / (X flops per loop)
    //
    // In this computation we ignore extra overhead that will be
    // an issue when local_dim is small.
    //

    double raw_epetra_time, thyra_wrapped_time;

    
    if(verbose) *out << "\n*** (C.1) Comparing the speed of RTOp verses raw Epetra_Vector operations\n";

    const double flop_adjust_factor_1 = 3.0;
    const int num_time_loops_1 = int( max_flop_rate / ( flop_adjust_factor_1 * local_dim * num_mv_cols ) ) + 1;

    {
        
      // Get references to Epetra_MultiVector objects in eV1 and eV2
      const RCP<Epetra_MultiVector>       eeV1 = get_Epetra_MultiVector(*epetra_map,eV1);
      const RCP<const Epetra_MultiVector> eeV2 = get_Epetra_MultiVector(*epetra_map,eV2);
      
      if(verbose)
        *out << "\nPerforming eeV1 = eeV2 (using raw Epetra_MultiVector::operator=(...)) " << num_time_loops_1 << " times ...\n";
      timer.start(true);
      for(int k = 0; k < num_time_loops_1; ++k ) {
        *eeV1 = *eeV2;
      }
      timer.stop();
      raw_epetra_time = timer.totalElapsedTime();
      if(verbose) *out << "  total time = " << raw_epetra_time << " sec\n";

      // When this block ends and eeV1 goes *out of scope then eV1 is guaranteed to be updated!
    }
    
    if(verbose)
      *out << "\nPerforming eV1 = eV2 (using Thyra::SpmdMultiVectorBase::applyOp(...)) " << num_time_loops_1 << " times ...\n";
    timer.start(true);
    for(int k = 0; k < num_time_loops_1; ++k ) {
      Thyra::assign( eV1.ptr(), *eV2 );
    }
    timer.stop();
    thyra_wrapped_time = timer.totalElapsedTime();
    if(verbose) *out << "  total time = " << thyra_wrapped_time << " sec\n";
    
    print_performance_stats( num_time_loops_1, raw_epetra_time, thyra_wrapped_time, verbose, *out );

    // RAB: 2004/01/05: Note, the above relative performance is likely
    // to be the worst of all of the others since RTOp operators are
    // applied seperately column by column but the relative
    // performance should go to about 1.0 when local_dim is
    // sufficiently large!  However, because
    // Epetra_MultiVector::Thyra::Assign(...) is implemented using double
    // pointer indexing, the RTOp implementation used with the Thyra
    // adapters is actually faster in some cases.  However, the extra
    // overhead of RTOp is much worse for very very small (order 10)
    // sizes.


    if(verbose)
      *out
        << "\n*** (C.2) Comparing Thyra::SpmdMultiVectorBase::apply() verses raw Epetra_MultiVector::Multiply()\n";

    Teuchos::TimeMonitor::zeroOutTimers();

    const double flop_adjust_factor_2 = 2.0;
    const int num_time_loops_2 = int( max_flop_rate / ( flop_adjust_factor_2* local_dim * num_mv_cols * num_mv_cols ) ) + 1;

    {
      
      // Get constant references to Epetra_MultiVector objects in eV1 and eV2
      const RCP<const Epetra_MultiVector> eeV1 = get_Epetra_MultiVector(*epetra_map,eV1);
      const RCP<const Epetra_MultiVector> eeV2 = get_Epetra_MultiVector(*epetra_map,eV2);
      
      Epetra_LocalMap eT_map((int) T->range()->dim(),0,*epetra_comm);
      Epetra_MultiVector eT(eT_map,T->domain()->dim());
      
      if(verbose)
        *out << "\nPerforming eeV1'*eeV2 (using raw Epetra_MultiVector::Multiply(...)) "	<< num_time_loops_2 << " times ...\n";
      timer.start(true);
      for(int k = 0; k < num_time_loops_2; ++k ) {
        eT.Multiply( 'T', 'N', 1.0, *eeV1, *eeV2, 0.0 );
      }
      timer.stop();
      raw_epetra_time = timer.totalElapsedTime();
      if(verbose) *out << "  total time = " << raw_epetra_time << " sec\n";
      
    }
    
    if(verbose)
      *out << "\nPerforming eV1'*eV2 (using Thyra::SpmdMultiVectorBase::apply(...)) "	<< num_time_loops_2 << " times ...\n";
    timer.start(true);
    for(int k = 0; k < num_time_loops_2; ++k ) {
      apply( *eV1, TRANS, *eV2, T.ptr() );
    }
    timer.stop();
    thyra_wrapped_time = timer.totalElapsedTime();
    if(verbose) *out << "  total time = " << thyra_wrapped_time << " sec\n";
  
    print_performance_stats( num_time_loops_2, raw_epetra_time, thyra_wrapped_time, verbose, *out );
    
    // RAB: 2004/01/05: Note, even though the Thyra adapter does
    // not actually call Epetra_MultiVector::Multiply(...), the
    // implementation in Thyra::SpmdMultiVectorBase::apply(...)
    // performs almost exactly the same flops and calls dgemm(...)
    // as well.  Herefore, except for some small overhead, the raw
    // Epetra and the Thyra wrapped computations should give
    // almost identical times in almost all cases.


    if(verbose) *out << "\n*** (C.3) Comparing Thyra::EpetraLinearOp::apply() verses raw Epetra_Operator::apply()\n";

    Teuchos::TimeMonitor::zeroOutTimers();

    const double flop_adjust_factor_3 = 10.0; // lots of indirect addressing
    const int num_time_loops_3 = int( max_flop_rate / ( flop_adjust_factor_3 * local_dim * num_mv_cols ) ) + 1;

    {
      
      // Get constant references to Epetra_MultiVector objects in eV1 and eV2
      const RCP<const Epetra_MultiVector> eeV1 = get_Epetra_MultiVector(*epetra_map,eV1);
      const RCP<Epetra_MultiVector>       eeY  = get_Epetra_MultiVector(*epetra_map,eY);
      
      if(verbose)
        *out << "\nPerforming eeY = eOp*eeV1 (using raw Epetra_Operator::apply(...)) " << num_time_loops_3 << " times ...\n";

      Teuchos::TimeMonitor::zeroOutTimers();

      timer.start(true);
      epetra_op->SetUseTranspose(false);
      for(int k = 0; k < num_time_loops_3; ++k ) {
        epetra_op->Apply( *eeV1, *eeY );
        //eeY->Scale(2.0);
      }
      timer.stop();

      raw_epetra_time = timer.totalElapsedTime();
      if(verbose) *out << "  total time = " << raw_epetra_time << " sec\n";

      if(verbose)
        Teuchos::TimeMonitor::summarize(*out << "\n");
      
    }
    
    if(verbose)
      *out << "\nPerforming eY = Op*eV1 (using Thyra::EpetraLinearOp::apply(...)) " << num_time_loops_3 << " times ...\n";

    Teuchos::TimeMonitor::zeroOutTimers();

    timer.start(true);
    for(int k = 0; k < num_time_loops_3; ++k ) {
      apply( *Op, NOTRANS, *eV1, eY.ptr() );
    }
    timer.stop();

    thyra_wrapped_time = timer.totalElapsedTime();
    if(verbose) *out << "  total time = " << thyra_wrapped_time << " sec\n";

    if(verbose)
      Teuchos::TimeMonitor::summarize(*out << "\n");
    
    print_performance_stats( num_time_loops_3, raw_epetra_time, thyra_wrapped_time, verbose, *out );

    // RAB: 2004/01/05: Note, the above Epetra adapter is a true
    // adapter and simply calls Epetra_Operator::apply(...) so, except
    // for some small overhead, the raw Epetra and the Thyra wrapped
    // computations should give ab*out exactly the same runtime for
    // almost all cases.

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if(verbose) {
    if(success) *out << "\nCongratulations! All of the tests seem to have run sucessfully!\n";
    else        *out << "\nOh no! at least one of the tests did not check out!\n";
  }

  return (success ? 0 : -1);

} // end main()
