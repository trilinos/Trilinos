// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// 

#include "Belos_LinearProblem.hpp"
#include "Belos_SummaryOutputterStatusTest.hpp"
#include "Belos_NativeNormStatusTest.hpp"
#include "Belos_ResidualNormStatusTest.hpp"
#include "TSFCoreDiagonalLinearOp.hpp"
#include "TSFCoreExplicitVectorView.hpp"
//#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

///
/** Class for generating a linear problem with RHS and
 * preconditioners, tolerances and stataus test.
 */
template<class Scalar>
class SimpleLinearProblemGenerator {
private:

	// ////////////////////////
	// Private types

	typedef Teuchos::ScalarTraits<Scalar>  ST;
	typedef typename ST::magnitudeType ScalarMag;

	enum { NUM_PREC_TYPE_OPTIONS = 4 };
	enum EPrecType { PREC_NONE, PREC_LEFT, PREC_RIGHT, PREC_LEFT_RIGHT };

	// ////////////////////////
	// Private data members

	int           dim_;
	double        opCondNum_;
	EPrecType     precType_;
	double        combinedOpCondNum_;
	int           numRhs_;
	int           blockSize_;
	ScalarMag     baseTol_;
	ScalarMag     varTolMag_;
	bool          useNativeNormStatusTest_;

public:

	// ////////////////////////
	// Public members

	/// Construct
	SimpleLinearProblemGenerator()
		:dim_(16),opCondNum_(10.0),precType_(PREC_NONE),combinedOpCondNum_(2.0),numRhs_(1),blockSize_(1)
		,baseTol_(1e-5),varTolMag_(0.0),useNativeNormStatusTest_(true)
		{}

	/// Setup CommandLineProcessor
	void setupCLP(
		bool                            *verbose
		,bool                           *dumpAll
		,Teuchos::CommandLineProcessor  *clp
		)
		{
			TEST_FOR_EXCEPT(clp==NULL);
			clp->setOption( "verbose", "quiet", verbose, "Determines if any output is printed or not" );
			clp->setOption( "dump-all", "no-dump", dumpAll, "Determines if quantities are dumped or not" );
			clp->setOption( "dim", &dim_, "Dimension of the linear system" );
			clp->setOption( "op-cond-num", &opCondNum_, "Condition number of the diagonal operator" );
			const EPrecType precTypeValues[NUM_PREC_TYPE_OPTIONS] = { PREC_NONE, PREC_LEFT, PREC_RIGHT, PREC_LEFT_RIGHT };
			const char*     precTypeNames[NUM_PREC_TYPE_OPTIONS]  = { "none",    "left",    "right", "left-right" };
			clp->setOption(
				"prec-type", &precType_
				,NUM_PREC_TYPE_OPTIONS, precTypeValues, precTypeNames
				,"Determines what type of preconditioner to use"
				);
			clp->setOption( "combined-op-cond-num", &combinedOpCondNum_, "Condition number of the preconditioned diagonal operator (if prec is used)" );
			clp->setOption( "num-rhs", &numRhs_, "Number of right-hand sides" );
			clp->setOption( "block-size", &blockSize_, "Number of RHSs solve simultaneously at a time" );
			clp->setOption( "base-tol", &baseTol_, "Base tolerance for random tolerances." );
			clp->setOption( "var-tol-mag", &varTolMag_, "Variation magnitude from base tolerance used for random tolerances." );
			clp->setOption( "use-native-status-test", "use-non-native-status-test", &useNativeNormStatusTest_
											, "Use the native or non-native residual status tests" );
		}
		
	/// Create tolerances and linear problem
	void setupLinearProblemEtc(
		const Teuchos::RefCountPtr<std::ostream>  &out
		,const bool                               verbose
		,const bool                               dumpAll
		,std::vector<ScalarMag>                   *tols
		,Belos::LinearProblemSetup<Scalar>        *lpup
		) const
		{
			TEST_FOR_EXCEPT(lpup==NULL||tols==NULL);
			using Teuchos::RefCountPtr;
			using Teuchos::rcp;
			//
			// A) Create operator, (optional) preconditioner, RHS and LHS
			//
			if(verbose)
				*out << "\nCreating diagonal opeator of dimension " << dim_ << " with condition number " << opCondNum_ << " ...\n";
			RefCountPtr<TSFCore::VectorSpace<Scalar> > space = rcp(new TSFCore::SerialVectorSpace<Scalar>(dim_));
			RefCountPtr<TSFCore::Vector<Scalar> > A_diag = space->createMember();
			if(1) {
				TSFCore::ExplicitMutableVectorView<Scalar> A_diag_ev(*A_diag);
				for( int k = 1; k <= dim_; ++k ) A_diag_ev(k) = ( opCondNum_ - 1 ) / ( dim_ - 1 ) * (k - 1) + 1;
			}
	    RefCountPtr<TSFCore::LinearOp<Scalar> > A = rcp( new TSFCore::DiagonalLinearOp<Scalar>(A_diag) );
			RefCountPtr<TSFCore::LinearOp<Scalar> > P_L, P_R;
			if( precType_!=PREC_NONE ) {
				double precCondNum = opCondNum_ / combinedOpCondNum_; 
				if( precType_==PREC_LEFT_RIGHT ) {
					precCondNum *= 0.5;
					if(verbose)
						*out
							<< "\nCreating diagonal left and right preconditioners of dimension "
							<< dim_ << " each with condition numbers " << precCondNum << " ...\n";
				}
				else {
					if(verbose)
						*out
							<< "\nCreating diagonal " << ( precType_==PREC_LEFT ? "left" : "right" )
							<< " preconditioner of dimension  " << dim_ << " with condition number " << precCondNum << " ...\n";
				}
				RefCountPtr<TSFCore::Vector<Scalar> > P_diag = space->createMember();
				if(1) {
					TSFCore::ExplicitMutableVectorView<Scalar> P_diag_ev(*P_diag);
					for( int k = 1; k <= dim_; ++k ) P_diag_ev(k) = 1 / ( ( precCondNum - 1 ) / ( dim_ - 1 ) * (k - 1) + 1 );
				}
				if( precType_==PREC_LEFT_RIGHT ) {
					P_L = rcp( new TSFCore::DiagonalLinearOp<Scalar>(P_diag) );
					P_R = rcp( new TSFCore::DiagonalLinearOp<Scalar>(P_diag) );
				}
				else {
					( precType_==PREC_LEFT ? P_L : P_R ) = rcp( new TSFCore::DiagonalLinearOp<Scalar>(P_diag) );
				}
				if(verbose)
					*out
						<< "\nCondition number of combined preconditioned operator is " << combinedOpCondNum_ << "!\n";
			}
			if(verbose)
				*out << "\nCreating RHS and LHS multi-vectors with " << dim_ << " rows and " << numRhs_ << " columns ...\n";
			RefCountPtr<TSFCore::MultiVector<Scalar> > B = space->createMembers(numRhs_),	X = space->createMembers(numRhs_);
			TSFCore::randomize( -ST::one(), +ST::one(), &*B );
			Teuchos::set_extra_data( space, "space", &X );
			Teuchos::set_extra_data( space, "space", &B );
			//
			// B) Create tolerances
			//
			tols->resize(numRhs_);
			ST::seedrandom(0);
			for( int k = 0; k < numRhs_; ++k ) (*tols)[k] = ST::magnitude( baseTol_ * std::pow(10,varTolMag_*ST::random()));
			if(verbose) {
				*out << "\nRandom tolerances: (j,tol) = {";
				for( int k = 0; k < numRhs_; ++k ) {
					if(k!=0) *out << ",";
					*out << "(" << k+1 << "," << (*tols)[k] << ")";
				}
				*out << "}\n";
			}
			//
			// C) Setup status test
			//
			RefCountPtr<Belos::AttachStatusTestBase<Scalar> > statusTest;
			// C.1) Summary outputting
			if(verbose) {
				statusTest = rcp(
					new Belos::SummaryOutputterStatusTest<Scalar>(
						out,std::string("  "),false
						,useNativeNormStatusTest_,!useNativeNormStatusTest_
						)
					);
			}
			// C.2) Native norm based
			RefCountPtr<Belos::AttachStatusTestBase<Scalar> > compStatusTest;
			if(useNativeNormStatusTest_)
				compStatusTest = rcp(new Belos::NativeNormStatusTest<Scalar>(numRhs_,&(*tols)[0]));
			else
				compStatusTest = rcp(new Belos::ResidualNormStatusTest<Scalar>(numRhs_,&(*tols)[0]));
			if(statusTest.get())
				statusTest->attach(compStatusTest);
			else
				statusTest = compStatusTest;
			// D) Setup the LinearProblem object
			lpup->setOperator(TSFCore::LinOpPersisting<Scalar>(A),Belos::OP_SYMMETRIC);
			if(P_L.get()) lpup->setLeftPrec(TSFCore::LinOpPersisting<Scalar>(P_L),Belos::OP_SYMMETRIC);
			if(P_R.get()) lpup->setRightPrec(TSFCore::LinOpPersisting<Scalar>(P_R),Belos::OP_SYMMETRIC);
			lpup->setRhs(B);
			lpup->setLhs(X);
			lpup->setBlockSize(blockSize_);
			lpup->setStatusTest(statusTest);
			lpup->completeSetup();
		}
};

///
/** Check the residual.
 */
template<class Scalar>
bool checkResidual(
	const TSFCore::LinOpNonPersisting<Scalar>                                  &A
	,const TSFCore::MultiVector<Scalar>                                        &B
	,const TSFCore::MultiVector<Scalar>                                        &X
	,const bool                                                                verbose
	,const std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>  &tols
	,std::ostream                                                              &out
	)
{
	typedef Teuchos::ScalarTraits<Scalar>  ST;
	typedef typename ST::magnitudeType ScalarMag;
	using Teuchos::CommandLineProcessor;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
	if(verbose)
		out << "\nChecking actual residual norms ...\n";
	const int numRhs = B.domain()->dim();
	RefCountPtr<TSFCore::MultiVector<Scalar> >
		R = B.range()->createMembers(numRhs);
	assign( &*R, B );
	A.apply( TSFCore::NOTRANS, X, &*R, ST::one(), -ST::one() );
	std::vector<ScalarMag> R_norms(numRhs);
	std::vector<ScalarMag> B_norms(numRhs);
	TSFCore::norms( *R, &R_norms[0] );
	TSFCore::norms( B, &B_norms[0] );
	bool success = true;
	for( int k = 0; k < numRhs; ++k ) {
		const ScalarMag R_norm_rel = R_norms[k] / B_norms[k];
		const bool result = ( R_norm_rel <= tols[k] );
		if(!result) success = false;
		if(verbose) {
			out
				<< "  ||R(:,"<<k+1<<")|| / ||B(:,j)|| = "<<R_norms[k]<<" / "<<B_norms[k]<<" = "<<R_norm_rel
				<< ( result ? " <= " : " > " ) << tols[k] << " : "
				<< ( result ? "passed" : "failed" )
				<< std::endl;
		}
	}
	return success;
}
