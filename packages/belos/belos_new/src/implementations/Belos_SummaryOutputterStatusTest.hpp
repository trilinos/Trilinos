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

#ifndef BELOS_SUMMARY_OUTPUTTER_STATUS_TEST_HPP
#define BELOS_SUMMARY_OUTPUTTER_STATUS_TEST_HPP

#include "Belos_AttachStatusTestBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "StandardCompositionMacros.hpp"

namespace Belos {

///
/** A status test subclass that print iteration summary
 * information only.
 *
 * This class will always just applies and returns the status from any
 * attached status test object (if a client called <tt>attach()</tt>).
 *
 * This class prints a running summary of each iteration, the maximum
 * native norm and (optionally) the native norm for each RHS.
 */
template<class Scalar>
class SummaryOutputterStatusTest : public AttachStatusTestBase<Scalar> {
public:

	/// Set the stream that output will be sent to
	STANDARD_COMPOSITION_MEMBERS( std::ostream, out );

	/// Set the leading string that will be printed at the beginning of each new line of output.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, leadingOutputStr );

	/// Print max norm over all current systems
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, printMaxNativeRhsNorm );

	/// Print norms for each native RHS separately.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, printEachNativeRhsNorm );

	/// Print norms for each original RHS separately.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, printEachOrigRhsNorm );

	/// Print norms for each original RHS numerator separately.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, printEachOrigRhsNormNumer );

	/// Print norms for each original RHS denominator separately.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, printEachOrigRhsNormDecom );

	/** @name Constructors / initializers */
	//@{
	
	///
	/** Construct to uninitialized.
	 *
	 * Postcondition:<ul>
	 * <li><tt>this->attachmentMode() == ATTACHED_TEST_INSTEAD</tt>
	 * <li><tt>this->get_out().get() == NULL</tt>
	 * <li><tt>this->leadingOutputStr() == ""</tt>
	 * <li><tt>this->printMaxNativeRhsNorm() == true</tt>
	 * <li><tt>this->printEachNativeRhsNorm() == false</tt>
	 * <li><tt>this->printEachOrigRhsNorm() == false</tt>
	 * </ul>
	 */
	SummaryOutputterStatusTest();

	///
	/** Construct with output stream and leading string optionally set.
	 *
	 * @param  out    [in] Smart pointer to output stream.
	 * @param  leadingOutputStr
	 *                [in] String that is printed at beginning of each new line.
	 *                Default <tt>leadingOutputStr==""</tt>.
	 * @param  printMaxNativeRhsNorm
	 *                [in] Determines if the max over all current systems is printed
	 *                or not.
	 *                Default <tt>printMaxNativeRhsNorm==true</tt>.
	 * @param  printEachNativeRhsNorm
	 *                [in] Determines if each native residual norm is printed as well or not.
	 *                Default <tt>printEachNativeRhsNorm==false</tt>.
	 * @param  printEachOrigRhsNorm
	 *                [in] Determines if each original residual norm is printed as well or not.
	 *                Default <tt>printEachOrigRhsNorm==false</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt></tt>
	 * </ul>
	 *
	 * Postcondition:<ul>
	 * <li><tt>this->attachmentMode() == ATTACHED_TEST_INSTEAD</tt>
	 * <li><tt>this->get_out().get() == out.get()</tt>
	 * <li><tt>this->leadingOutputStr() == leadingOutputStr</tt>
	 * <li><tt>this->printMaxNativeRhsNorm() == printMaxNativeRhsNorm</tt>
	 * <li><tt>this->printEachNativeRhsNorm() == printEachNativeRhsNorm</tt>
	 * <li><tt>this->printEachOrigRhsNorm() == printEachOrigRhsNorm</tt>
	 * </ul>
	 */
	SummaryOutputterStatusTest(
		const Teuchos::RefCountPtr<std::ostream> &out
		,const std::string                       &leadingOutputStr       = std::string("")
		,const bool                              printMaxNativeRhsNorm   = true
		,const bool                              printEachNativeRhsNorm  = false
		,const bool                              printEachOrigRhsNorm    = false
		,const bool                              printEachOrigRhsNormNumer = false
		,const bool                              printEachOrigRhsNormDecom = false
		);

	//@}

	/** @name Overridden from StatusTest */
	//@{

	///
	void reset();
	///
	void checkStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		);

	//@}

protected:

	/** @name Overridden from AttachedStatusTestBase */
	//@{

	/// Never called
	void protectedReset();
	/// Never called
	void protectedCheckStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		);

	//@}

private:

	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMag;

	bool resetCalled_;
	std::valarray<ScalarMag>  R_native_norms_;
	std::valarray<ScalarMag>  R_bar_norms_;
	std::valarray<ScalarMag>  R_bar_0_norms_;
	std::valarray<int>        currRhsIndexes_;
	int                       currNumRestarts_;

};

// //////////////////////////////
// Implementations


// Constructors / initializers

template<class Scalar>
SummaryOutputterStatusTest<Scalar>::SummaryOutputterStatusTest()
	:AttachStatusTestBase<Scalar>(ATTACHED_TEST_INSTEAD)
	,printMaxNativeRhsNorm_(true)
	,printEachNativeRhsNorm_(false)
	,printEachOrigRhsNorm_(false)
	,printEachOrigRhsNormNumer_(false)
	,printEachOrigRhsNormDecom_(false)
	,resetCalled_(true)
	,currNumRestarts_(0)
{}

template<class Scalar>
SummaryOutputterStatusTest<Scalar>::SummaryOutputterStatusTest(
	const Teuchos::RefCountPtr<std::ostream>  &out
	,const std::string                        &leadingOutputStr
	,const bool                               printMaxNativeRhsNorm
	,const bool                               printEachNativeRhsNorm
	,const bool                               printEachOrigRhsNorm
	,const bool                               printEachOrigRhsNormNumer
	,const bool                               printEachOrigRhsNormDecom
	)
	:AttachStatusTestBase<Scalar>(ATTACHED_TEST_INSTEAD)
	,out_(out)
	,leadingOutputStr_(leadingOutputStr)
	,printMaxNativeRhsNorm_(printMaxNativeRhsNorm)
	,printEachNativeRhsNorm_(printEachNativeRhsNorm)
	,printEachOrigRhsNorm_(printEachOrigRhsNorm)
	,printEachOrigRhsNormNumer_(printEachOrigRhsNormNumer)
	,printEachOrigRhsNormDecom_(printEachOrigRhsNormDecom)
	,resetCalled_(true)
	,currNumRestarts_(0)
{
	TEST_FOR_EXCEPT(out.get()==NULL);
}

// Overridden from StatusTest

template<class Scalar>
void SummaryOutputterStatusTest<Scalar>::reset()
{
	Teuchos::RefCountPtr<StatusTest<Scalar> > attachedStatusTest = getAttachedStatusTest();
	if(attachedStatusTest.get()) attachedStatusTest->reset();
	resetCalled_ = true;
}

template<class Scalar>
void SummaryOutputterStatusTest<Scalar>::checkStatus(
	const BasicIterationState<Scalar>         &bis
	,const int                                currBlockSize
	,const int                                currNumRhs
	,EStatusType                              status[]
	)
{
	Teuchos::RefCountPtr<StatusTest<Scalar> > attachedStatusTest = getAttachedStatusTest();
	if(attachedStatusTest.get()) attachedStatusTest->checkStatus(bis,currBlockSize,currNumRhs,status);
	const LinearProblemState<Scalar> &lps = bis.getProblem();
	const bool getCurrRhsIndexes = resetCalled_ || printEachNativeRhsNorm() || printEachOrigRhsNorm();
	if(getCurrRhsIndexes) {
		if( static_cast<int>(currRhsIndexes_.size()) < currBlockSize ) currRhsIndexes_.resize(currBlockSize);
		lps.getCurrRhsIndexes( currNumRhs, &currRhsIndexes_[0] );
	}
	if(out_.get()) {
		const std::string &leadstr = leadingOutputStr();
		// Print header
		if(resetCalled_) {
			*out_
				<< leadstr << "Solving RHSs = ["
				<< currRhsIndexes_[0] << ",...," << currRhsIndexes_[currNumRhs-1] << "]"
				<< "\n  , solver = \'" << typeid(bis).name() << "\'"
				<< "\n  , attachedStatusTest = \'";
			if(attachedStatusTest.get()) *out_ << typeid(*attachedStatusTest).name();
			else *out_ << "NULL";
			*out_	
				<< "\'"
				<< "\n  , currNativeResidualType = \'" << toString(bis.getCurrNativeResidualType()) << "\'\n"; 
		}
		// Compute/get norms
		if(printMaxNativeRhsNorm() || printEachNativeRhsNorm()) {
			if( static_cast<int>(R_native_norms_.size()) < currBlockSize ) R_native_norms_.resize(currBlockSize);
			bis.getCurrNativeResiduals(currBlockSize,&R_native_norms_[0]);
		}
		if(printEachOrigRhsNorm()) {
			if(!lps.isCurrLhsUpdated()) bis.forceCurrLhsUpdate();	
			const TSFCore::MultiVector<Scalar>
				&R_bar_0 = lps.getCurrInitResidual(),
				&R_bar   = lps.getCurrResidual();
			if(static_cast<int>(R_bar_norms_.size()) < currBlockSize) {
				R_bar_norms_.resize(currBlockSize);
				R_bar_0_norms_.resize(currBlockSize);
			}
			TSFCore::norms( R_bar, &R_bar_norms_[0] );
			TSFCore::norms( R_bar_0, &R_bar_0_norms_[0] );
		}
		// Print output for this iteration
		const int currNumRestarts = bis.getCurrNumRestarts();
		if( currNumRestarts_ != currNumRestarts ) {
			*out_<<leadstr<<"Restart performed (currNumRestarts="<<currNumRestarts<<")\n";
			currNumRestarts_ = currNumRestarts;
		}
		*out_ << leadstr << "  iter="<<bis.getCurrNumIters()<<", currNumRhs="<<currNumRhs;
		if(printMaxNativeRhsNorm()) {
			const Scalar maxNorm = *std::max_element(&R_native_norms_[0],&R_native_norms_[0]+currBlockSize);
			*out_<<"max{||R_native||/||R_native_0||}="<<maxNorm<<std::endl;
		}
		else {
			*out_<<std::endl;
		}
		if( printEachNativeRhsNorm() ) {
			*out_ << leadstr << "    (j(k),||R_native(:,k)||/||R_native_0(:,k)||) = {";
			for( int k = 0; k < currNumRhs; ++k ) {
				if( k != 0 ) *out_ << ",";
				*out_ << "(" << currRhsIndexes_[k] << "," << R_native_norms_[k] << ")";
			}
			*out_ << "}\n";
		}
		if( printEachOrigRhsNorm() ) {
			*out_ << leadstr << "    (j(k),||R_orig(:,k)||/||R_orig_0(:,k)||) = {";
			for( int k = 0; k < currNumRhs; ++k ) {
				if( k != 0 ) *out_ << ",";
				*out_ << "(" << currRhsIndexes_[k] << "," << (R_bar_norms_[k]/R_bar_0_norms_[k]) << ")";
			}
			*out_ << "}\n";
		}
		if( printEachOrigRhsNormNumer() ) {
			*out_ << leadstr << "    (j(k),||R_orig(:,k)||) = {";
			for( int k = 0; k < currNumRhs; ++k ) {
				if( k != 0 ) *out_ << ",";
				*out_ << "(" << currRhsIndexes_[k] << "," << (R_bar_norms_[k]) << ")";
			}
			*out_ << "}\n";
		}
		if( printEachOrigRhsNormDecom() ) {
			*out_ << leadstr << "    (j(k),||R_orig_0(:,k)||) = {";
			for( int k = 0; k < currNumRhs; ++k ) {
				if( k != 0 ) *out_ << ",";
				*out_ << "(" << currRhsIndexes_[k] << "," << (R_bar_0_norms_[k]) << ")";
			}
			*out_ << "}\n";
		}
	}
	resetCalled_ = false;
}

// Overridden from AttachedStatusTestBase

template<class Scalar>
void SummaryOutputterStatusTest<Scalar>::protectedReset()
{
	TEST_FOR_EXCEPT(true); // Never called but just in case ...
}

template<class Scalar>
void SummaryOutputterStatusTest<Scalar>::protectedCheckStatus(
	const BasicIterationState<Scalar>         &bis
	,const int                                currBlockSize
	,const int                                currNumRhs
	,EStatusType                              status[]
	)
{
	TEST_FOR_EXCEPT(true); // Never called but just in case ...
}

} // namespace Belos

#endif  // BELOS_SUMMARY_OUTPUTTER_STATUS_TEST_HPP
