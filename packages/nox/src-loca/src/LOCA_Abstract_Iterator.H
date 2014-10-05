// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_ABSTRACT_ITERATOR_H
#define LOCA_ABSTRACT_ITERATOR_H

namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  namespace Abstract {

    /*!
      \brief An abstract interface for implementing iteration.
    */
    /*!
      The %LOCA::Abstract::Iterator defines an interface for implementing
      many kinds of iterative processes.  In %LOCA, this is used to implement
      the %Stepper which computes points along a continuation curve.

      Many iterative processes can be abstracted in the following manner:

      <ul>
      <li> Initialize iteration (start)
      <li> Compute iteration (iterate)
           <ul>
           <li> while iterator is not finished
                <ul>
                <li> preprocess step (preprocess)
                <li> compute step (compute)
            <li> postprocess step (posprocess)
            <li> check iterator status (stop)
            </ul>
            </ul>
      <li> Finalize iteration (finish)
      </ul>

      The run method of the iterator implements this iterative process with
      start, finish, preprocess, compute and postprocess left as pure
      virtual methods to be implemented for the specific iterative process.

      The iterator has one parameter, "Max Steps" (default 100) giving the
      maximum number of steps the iterator should take.  The default
      implementation of stop only stops the iterator when this maximum
      number of steps has been reached.
    */
    class Iterator {

    public:

      //! Enumerated type for status of the iterator
      enum IteratorStatus {
        FinishedWithNatural = 4,  //!< A successful last step was taken to hit the target
        BeforeFirstIteration = 3, //!< The iteration has not started
        LastIteration = 2,        //!< This is the last iteration
        Finished = 1,             //!< The iterator is finished
        Failed = 0,               //!< The iterator has failed
        NotFinished = -1          //!< The iterator is not finished
      };

      //! Enumerated type for status of each step of iterator
      enum StepStatus {
        Successful = 1,    //!< The last step was successful
        Unsuccessful = 0,  //!< The last step was unsuccessful
        Provisional = -1   //!< One more unsuccessful step means failure
      };

      //! Constructor
      Iterator(Teuchos::ParameterList& p);

      //! Copy Constructor
      Iterator(const Iterator& it);

      //! Destructor
      virtual ~Iterator();

      //! Reset the iterator to start a new iteration
      virtual bool resetIterator(Teuchos::ParameterList& p);

      //! Return the status of the iterator
      virtual IteratorStatus getIteratorStatus() const;

      //! Returns the number of accepted steps
      virtual int getStepNumber() const;

      //! Returns the number of failed steps
      virtual int getNumFailedSteps() const;

      //! Returns the total number of steps attempted
      virtual int getNumTotalSteps() const;

      //! Run the iterator
      virtual IteratorStatus run();

      /// \brief Initialize iterating one by one.
      ///
      /// \param iteratorStatus [out] On exit, the value to be passed to \c
      /// iterateOnce().
      ///
      /// Instead of running all iterations at once using \c run(), iterate one
      /// by one. First initialize the iterator using this method.
      virtual void initializeIterateOnce(IteratorStatus& iteratorStatus);

      /// \brief Iterate once.
      ///
      /// \param iteratorStatus [in, out] On entry, either the value from
      /// initializeIterateOnce() or from the previous call of iterateOnce(), if
      /// it returned true. On exit, the status of the iterator. If
      /// iterateOnce() returns true, also the input for the next call.
      ///
      /// \param stepStatus [in, out] Same behavior as iteratorStatus, but for
      /// StepStatus.
      ///
      /// Returns \c true if one should keep iterating, false if done.
      ///
      /// Typical usage is as follows:
      /// \code
      /// Iterator::IteratorStatus istat;
      /// Iterator::StepStatus sstat;
      /// it->initializeIterateOnce(istat);
      /// while (it->iterate(istat, sstat)) {
      ///   // Do something with the data now that we are at the end of
      ///   // this continuation step.
      /// }
      /// // Do something different with the data at the last continuation step.
      /// \endcode
      /// Alternatively, if you do not want to distinguish between intermediate
      /// and the final continuation steps, you can structure the loop this way:
      /// \code
      /// it->initializeIterateOnce(istat);
      /// for (;;) {
      ///   const bool keep_iterating = it->iterate(istat, sstat);
      ///   // Do something with the data now that we are at the end of
      ///   // this continuation step.
      ///   if (!keep_iterating) break;
      /// }
      /// \endcode
      virtual bool iterateOnce(
        IteratorStatus& iteratorStatus, StepStatus& stepStatus);

    protected:

      //! Constructor
      Iterator();

      //! Perform iteration
      virtual IteratorStatus iterate();

      //! Check stopping criteria
      virtual IteratorStatus stop(StepStatus);

      //! Intialize iterator
      virtual IteratorStatus start() = 0;

      //! Finalize iterator
      virtual IteratorStatus finish(IteratorStatus) = 0;

      //! Preprocess step
      virtual StepStatus preprocess(StepStatus) = 0;

      //! Compute step
      virtual StepStatus compute(StepStatus) = 0;

      //! Postprocess step
      virtual StepStatus postprocess(StepStatus) = 0;

      //! Set iterator status as LastIteration
      virtual void setLastIteration();

      //! Returns true of the iteratorStatus is LastIteration
      virtual bool isLastIteration();

      //! Compute step, iterator status from status of each component of step
      virtual StepStatus computeStepStatus(StepStatus preStatus,
                       StepStatus compStatus,
                       StepStatus postStatus);

    protected:

      //! Current step number of continuation algorithm (does not include failed steps).
      int stepNumber;

      //! Number of continuation steps that have failed.
      int numFailedSteps;

      //! Total number of steps attempeted (includes failed and successful steps).
      int numTotalSteps;

      //! Maximum number of continuation steps to take. Defaults to 100.
      int maxSteps;

      //! Current status of the iterator
      IteratorStatus iteratorStatus;

    };

  } // end namespace Abstract
} // end namespace LOCA

#endif

