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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_UTILS_COUNTERS_H
#define NOX_LINESEARCH_UTILS_COUNTERS_H

#include "NOX_Common.H"

// Forward declaration
namespace NOX {
  class Utils;
  namespace Parameter {
    class List;
  }
}

namespace NOX {

namespace LineSearch {

//! %Common line search utilites
namespace Utils {

//! Common counters that all line search algorithms should report.
/*!  

  <B> Output list keys </B> <br>

Line searches have a common set of data that should be tracked and saved in the parameter list for the users.  This class provides a set of common data objects, accessors, and a routine to print them.  A sublist for output parameters will be created called "Output" in the "Line Search" parameter sublist.  Valid output keys for the parameter list are:

  - "Total Number of Line Search Calls" - Total number of calls to the compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of steps that could not directly take a full step and meet the required "Convergence Criteria" (i.e. The line search had to reduce the step length using inner iteration calculations over iterate \f$ k \f$). 

  - "Total Number of Failed Line Searches" - total number of line searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of inner iterations \f$ k \f$ performed by this object.

*/

class Counters {

public:

  //! Default constructor.
  Counters();

  //! Destructor.
  virtual ~Counters();

  //! Reset the counters .
  virtual void reset();

  //! Sets the common line search data in an "Output" sublist of the lineSearchParams list that is supplied in the method call. 
  virtual bool setValues(NOX::Parameter::List& lineSearchParams);

  /** @name Increment Methods
      

  */ 
  //@{

  //! Increase the counter for the total number of line search calls by numLS.
  virtual void increaseNumLineSearches(int numLS=1);

  //! Increase the counter for the total number of non-trivial line search calls by numNonTrivialLS.
  virtual void increaseNumNonTrivialLineSearches(int numNonTrivialLS=1);

  //! Increase the counter for the total number of failed line search calls by numFailedLS
  virtual void increaseNumFailedLineSearches(int numFailedLS=1);

  //! Increase the counter for the total number of line search inner iterations performed by the line search.
  virtual void increaseNumIterations(int numInnerIters=1);

  //@}

  /** @name Accessor Methods
      
  Returns the current counter value 

  */ 
  //@{

  //! Increase the counter for the total number of line search calls by numLS.
  virtual int getNumLineSearches() const;

  //! Increase the counter for the total number of non-trivial line search calls by numNonTrivialLS.
  virtual int getNumNonTrivialLineSearches() const;

  //! Increase the counter for the total number of failed line search calls by numFailedLS
  virtual int getNumFailedLineSearches() const;

  //! Increase the counter for the total number of line search inner iterations performed by the line search.
  virtual int getNumIterations() const;

  //@}


private:

  /** @name Common Line Search Data Members 
  
  All line searches in general should report the following information.  We save a pointer to the parameter list associated with the line search and set an output sublist with each of the following parameters.
  */ 
  //@{

  //! Counter for the total number of times the compute() method was called for this object. 
  int totalNumLineSearchCalls;

  //! Counter for number of times the line search was actually needed.  This is the number of times the compute() method was called for this object where a step length of 1.0 could not be used.
  int totalNumNonTrivialLineSearches;

  //! Counter for the number of times the compute() method was called where the line search algorithm failed.
  int totalNumFailedLineSearches;

  //! Counter for the total number of line search inner iterations performed.
  int totalNumIterations;

  //@}

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
