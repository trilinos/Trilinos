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

#ifndef BELOS_ATTACH_STATUS_TEST_BASE_HPP
#define BELOS_ATTACH_STATUS_TEST_BASE_HPP

#include "Belos_StatusTest.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Belos {

///
/** Node base subclass for all <tt>StatusTest</tt> subclasses
 * that want to allow the attachment of status tests.
 *
 * This subclass sets up the machinary for setting up an attached
 * status test and considering it in testing.
 *
 * <b>Note to subclass developers:</b>
 *
 * The only functions that must be overridden are the protected
 * functions <tt>protectedReset()</tt> and
 * <tt>protectedCheckStatus()</tt>.
 */
template<class Scalar>
class AttachStatusTestBase : public virtual StatusTest<Scalar> {
public:

	/** @name Public types */
	//@{

	///
	enum EAttachmentMode {
		ATTACHED_TEST_EXCLUDE  ///< The attached status test will never even be called.
		,ATTACHED_TEST_IGNORE  ///< The attached status test will be called but the result ignored.
		,ATTACHED_TEST_INSTEAD ///< The attached status test will be used in place of this status test.
		,ATTACHED_TEST_AND     ///< The attached status test will be called and have its results and'ed.
		,ATTACHED_TEST_OR      ///< The attached status test will be called and have its results or'ed
	};

	//@}

	///
	typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMagnitude;

	/** @name Constructors / initializers */
	//@{

	/// Calls <tt>attachmentMode()</tt>
	AttachStatusTestBase(
		const EAttachmentMode attachmentMode = ATTACHED_TEST_AND
		);

	/// Set the attachment mdoe (see <tt>EAttachmentMode</tt>)
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EAttachmentMode, attachmentMode );

	//@}

	/** @name StatusTest attachment */
	//@{

	///
	/** Attach another status test object.
	 *
	 * @param  statusTest  [in] Smart pointer to another status test that
	 *                     can be considered in some way.  It is allowed for
	 *                     <tt>statusTest.get()==NULL</tt> in which case
	 *                     the current status test (if one is currently
	 *                     attached) will be unattached.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getAttachedStatusTest().get() == statusTest.get()</tt>
	 * </ul>
	 *
	 * The convergence tester <tt>statusTest</tt> being attached can be
	 * dealt any way that <tt>*this</tt> chooses.
	 */
	virtual void attach( const Teuchos::RefCountPtr<StatusTest<Scalar> > &statusTest );

	///
	/** Get a smart pointer to non-<tt>const</tt> attached status test.
	 */
	virtual Teuchos::RefCountPtr<StatusTest<Scalar> > getAttachedStatusTest();

	///
	/** Get a smart pointer to <tt>const</tt> attached status test.
	 * tester.
	 *
	 * The default implementation returns
	 \code

	   const_cast<StatusTest<Scalar>*>(this)->getAttachedStatusTest()

   \endcode
	 *
	 * No override of this function should be needed.
	 */
	virtual Teuchos::RefCountPtr<const StatusTest<Scalar> > getAttachedStatusTest() const;

	//@}

	/** @name Overridden from StatusTest */
	//@{

	///
	/** Overridden.
	 *
	 * Overridden to call <tt>protectedReset()</tt> and will call
	 * <tt>getAttachedStatusTest()->reset()<tt> if
	 * <tt>getAttachedStatusTest().get()!=NULL &&
	 * attachmentMode()!=ATTACH_TEST_EXCLUDE</tt>.
	 */
	void reset();
	///
	/** Overridden.
	 *
	 * Overridden to call <tt>protectedCheckStatus()</tt> and will call
	 * <tt>getAttachedStatusTest()->checkStatus()<tt> if
	 * <tt>getAttachedStatusTest().get()!=NULL &&
	 * attachmentMode()!=ATTACH_TEST_EXCLUDE</tt>.
	 *
	 * The way that that status test returned form
	 * <tt>getAttachedStatusTest()->checkStatus()<tt> is dealt with is
	 * determined by the value returned by <tt>attachmentMode()</tt>
	 * just before this function is called.  It should be obvious what
	 * the post conditions are based on the value of
	 * <tt>attachmentMode()</tt>.
	 */
	void checkStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		);

	//@}

protected:

	/** @name Protected virtual functions to be overridden */
	//@{

	///
	virtual void protectedReset() {}

	///
	virtual void protectedCheckStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		) = 0;

	//@}

private:

	Teuchos::RefCountPtr<StatusTest<Scalar> > attachedStatusTest_;

	static EStatusType andOrStatus( const EStatusType status1, const EStatusType status2 );

};

// /////////////////////////////////////
// Implementations

// Constructors / initializers

template<class Scalar>
AttachStatusTestBase<Scalar>::AttachStatusTestBase(
	const EAttachmentMode attachmentMode
	)
{
	this->attachmentMode(attachmentMode);
}

// StatusTest attachment

template<class Scalar>
void AttachStatusTestBase<Scalar>::attach(
	const Teuchos::RefCountPtr<StatusTest<Scalar> > &statusTest
	)
{
	attachedStatusTest_ = statusTest;
}

template<class Scalar>
Teuchos::RefCountPtr<StatusTest<Scalar> >
AttachStatusTestBase<Scalar>::getAttachedStatusTest()
{
	return attachedStatusTest_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const StatusTest<Scalar> >
AttachStatusTestBase<Scalar>::getAttachedStatusTest() const
{
	return const_cast<AttachStatusTestBase<Scalar>*>(this)->getAttachedStatusTest();
}

// Overridden from StatusTest

template<class Scalar>
void AttachStatusTestBase<Scalar>::reset()
{
	protectedReset();
	if( attachedStatusTest_.get() && attachmentMode()!=ATTACHED_TEST_EXCLUDE )
		attachedStatusTest_->reset();
}

template<class Scalar>
void AttachStatusTestBase<Scalar>::checkStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
	)
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	protectedCheckStatus(bis,currBlockSize,currNumRhs,status);

	if( attachedStatusTest_.get() && attachmentMode()!=ATTACHED_TEST_EXCLUDE ) {
		Workspace<EStatusType> status2(wss,currNumRhs,false);
		attachedStatusTest_->checkStatus(bis,currBlockSize,currNumRhs,&status2[0]);
		switch(attachmentMode()) {
			case ATTACHED_TEST_IGNORE: {
				break; // Just ignore what is in status2
			}
			case ATTACHED_TEST_INSTEAD: {
				for( int k = 0; k < currNumRhs; ++k ) status[k] = status2[k];
				break;
			}
			case ATTACHED_TEST_AND: {
				for( int k = 0; k < currNumRhs; ++k ) {
					if( status[k] == STATUS_CONVERGED && status2[k] == STATUS_CONVERGED )
						status[k] = STATUS_CONVERGED;
					else
						status[k] = andOrStatus(status[k],status2[k]);
				}
				break;
			}
			case ATTACHED_TEST_OR: {
				for( int k = 0; k < currNumRhs; ++k ) {
					if( status[k] == STATUS_CONVERGED || status2[k] == STATUS_CONVERGED )
						status[k] = STATUS_CONVERGED;
					else
						status[k] = andOrStatus(status[k],status2[k]);
				}
				break;
			}
			default:
				TEST_FOR_EXCEPT(true);
		}
	}
}

// private

template<class Scalar>
EStatusType AttachStatusTestBase<Scalar>::andOrStatus( const EStatusType status1, const EStatusType status2 )
{
	if( status1 == STATUS_CONVERGED && status2 == STATUS_CONVERGED )
		return STATUS_CONVERGED;
	else if ( status1 == STATUS_NAN || status2 == STATUS_NAN )
		return STATUS_NAN;
	else if ( status1 == STATUS_FAILED || status2 == STATUS_FAILED )
		return STATUS_FAILED;
	else if ( status1 == STATUS_UNCONVERGED || status2 == STATUS_UNCONVERGED )
		return STATUS_UNCONVERGED;
	else
		return STATUS_UNCHECKED;
}

} // namespace Belos

#endif  // BELOS_ATTACH_STATUS_TEST_BASE_HPP
