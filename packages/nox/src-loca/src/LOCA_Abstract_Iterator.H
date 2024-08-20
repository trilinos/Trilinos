// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

