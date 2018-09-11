// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef NOX_LINESEARCH_UTILS_COUNTERS_H
#define NOX_LINESEARCH_UTILS_COUNTERS_H

#include "NOX_Common.H"

// Forward declaration
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;
}

namespace NOX {

//! Common counters that all line search algorithms should report.
/*!

  <B> Output list keys </B> <br>

Line searches have a common set of data that should be tracked and
saved in the parameter list for the users.  This class provides a set
of common data objects, accessors, and a routine to print them.  A
sublist for output parameters will be created called "Output" in the
"Line Search" parameter sublist.  Valid output keys for the parameter
list are:

  - "Total Number of Line Search Calls" - Total number of calls to the
    compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of
    steps that could not directly take a full step and meet the
    required "Convergence Criteria" (i.e. The line search had to
    reduce the step length using inner iteration calculations over
    iterate \f$ k \f$).

  - "Total Number of Failed Line Searches" - total number of line
    searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of
    inner iterations \f$ k \f$ performed by this object.

*/

class LineSearchCounters {

public:

  //! Default constructor.
  LineSearchCounters();

  //! Destructor.
  virtual ~LineSearchCounters();

  //! Reset the counters .
  virtual void reset();

  //! Sets the common line search data in an "Output" sublist of the
  //! lineSearchParams list that is supplied in the method call.
  virtual bool setValues(Teuchos::ParameterList& lineSearchParams);

  /** @name Increment Methods


  */
  //@{

  //! Increment the counter for the total number of line search calls by n
  virtual void incrementNumLineSearches(int n = 1);

  //! Increment the counter for the total number of non-trivial line search calls by n
  virtual void incrementNumNonTrivialLineSearches(int n = 1);

  //! Increment the counter for the total number of failed line search calls by n
  virtual void incrementNumFailedLineSearches(int n = 1);

  //! Increment the counter for the total number of line search inner iterations by n
  virtual void incrementNumIterations(int n = 1);

  //@}

  /** @name Accessor Methods

  Returns the current counter value

  */
  //@{

  //! Return the counter for the total number of line search calls
  virtual int getNumLineSearches() const;

  //! Return the counter for the total number of non-trivial line search calls
  virtual int getNumNonTrivialLineSearches() const;

  //! Return the counter for the total number of failed line search calls
  virtual int getNumFailedLineSearches() const;

  //! Return the counter for the total number of line search inner iterations
  virtual int getNumIterations() const;

  //@}


private:

  /** @name Common Line Search Data Members

  All line searches in general should report the following
  information.  We save a pointer to the parameter list associated
  with the line search and set an output sublist with each of the
  following parameters.
  */
  //@{

  //! Counter for the total number of times the compute() method was
  //! called for this object.
  int totalNumLineSearchCalls;

  //! Counter for number of times the line search was actually needed.
  //! This is the number of times the compute() method was called for
  //! this object where a step length of 1.0 could not be used.
  int totalNumNonTrivialLineSearches;

  //! Counter for the number of times the compute() method was called
  //! where the line search algorithm failed.
  int totalNumFailedLineSearches;

  //! Counter for the total number of line search inner iterations performed.
  int totalNumIterations;

  //@}

};

} // namespace NOX

#endif
