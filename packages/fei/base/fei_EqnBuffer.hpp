#ifndef _fei_EqnBuffer_hpp_
#define _fei_EqnBuffer_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fwd.hpp"

#include <vector>

/**
A class for holding equation data, along with optional RHS coefficients for the
equations. 
*/

class EqnBuffer {
 public:
  /** Constructor */
  EqnBuffer();

  /** copy constructor */
  EqnBuffer(const EqnBuffer& src);

  /** Destructor */
  virtual ~EqnBuffer();

   /** Return a 'clone' of this object, including a copy of all internal data.
    */
   EqnBuffer* deepCopy();

   /** assignment operator */
   EqnBuffer& operator=(const EqnBuffer& src);

   /** Return the number of equations held in this object.*/
   int getNumEqns() {return(eqnNumbers_.size());};

   /** Return a list of the equation-numbers held in this object. */
   std::vector<int>& eqnNumbers() {return(eqnNumbers_);};

   /** Return a table (actually an array of pointers to fei::CSVecs) of the
       equations. Number-of-arrays is 'getNumEqns', length of i-th array
       is 'lengthsPtr()[i]'. */
   std::vector<fei::CSVec*>& eqns() {return(eqns_);};

   /** Return the number of right-hand-side coefficientss stored for each
       equation.*/
   int getNumRHSs() {return(numRHSs_);};

   /** Set the number of right-hand-side coefficients to be stored for each
       equation. This function internally allocates a table to hold the rhs
       coefficients. This table has number-of-rows == 'getNumEqns()', and
       number-of-columns == 'n'. Note that this function should not be called
       after rhs coefficients have been added, because it destroys and
       re-allocates the rhs-coef table.
       @param n
   */
   void setNumRHSs(int n);

   /** Set a right-hand-side coefficient. 
       @param eqnNumber Must correspond to one of the equation-numbers in the
       list 'eqnNumbers()'.
       @param rhsIndex Must be in the range [0 .. getNumRHSs()-1].
       @param value The coefficient.
   */
   int addRHS(int eqnNumber, int rhsIndex, double value, bool accumulate=true);

   /** Return the table (actually an array of pointers to arrays) of
       right-hand-side coefficients. Number-of-arrays == 'getNumEqns()',
       number-of-columns == 'getNumRHSs()'.
   */
   std::vector<std::vector<double>*>* rhsCoefsPtr() {return(&rhsCoefs_);};

   /** Return an offset into the 'eqnNumbers()' list, being the position at
       which equation-number 'eqn' is located.
       @param eqn
       @return offset 'eqn's location, or -1 if 'eqn' is not found.
   */
   int getEqnIndex(int eqn);

   /** Query whether 'eqn' is present in the table of column-indices.
       @param eqn
       @return offset Row of the table that contains 'eqn', or -1 if not found.
   */
   int isInIndices(int eqn);

   /** Add an equation, with associated coefficients and column-indices.
       @param eqnNumber 
       @param coefs
       @param indices
       @param len Number of coefficients and indices
       @param accumulate If this equation is already present, indicate whether
       incoming coefficients should replace existing ones, or be accumulated 
       into them. Note that this only matters if eqnNumber is already present
       AND if one or more of 'indices' are already present.
       @param create_indices_union defaults to false. If true, adds each set
       of indices to a union of all indices. This will be used if the
       'isInIndices' method is called later.
       @return error-code, 0 if successful
   */
   int addEqn(int eqnNumber, const double* coefs, const int* indices,
               int len, bool accumulate, bool create_indices_union=false);

   /** Add the equations in inputEqns. Upon completion, the set of equations in
       'this' EqnBuffer is the union of the equations that were already present,
       and the equations in inputEqns.
       @param inputEqns EqnBuffer containing equations to be added to this
       EqnBuffer.
       @param accumulate If an equation is already present, determines whether
       to add its coefficients to those already in place, or replace them.
   */
   int addEqns(EqnBuffer& inputEqns, bool accumulate);

   /** Given an equation number and a column-index in that equation, get the
       corresponding coefficient.
       @param eqnNumber Input Equation-number
       @param colIndex Input Column-index in equation-number.
       @param coef Output coefficient corresponding to colIndex in eqnNumber.
       @return error-code 0 if successful. If eqnNumber is not found, or if
       equation 'eqnNumber' doesn't contain 'colIndex', then -1 is returned and
       coef is not referenced.
   */
   int getCoef(int eqnNumber, int colIndex, double& coef);

   /** Remove the specified column-index from the specified equation-number.
       @return error-code 0 if successful. If eqnNumber is not found, -1 is
       returned. If equation eqnNumber doesn't contain 'colIndex', the remove
       request is considered successful and 0 is returned.
   */
   int removeIndex(int eqnNumber, int colIndex);

   /** Combine the functions 'getCoef()' and 'removeIndex()'. Has the same
    effect as calling those two functions separately, but saves two searches.
    @return error-code 0 if successful. -1 is returned if position (eqnNumber,
    colIndex is not found.
   */
   int getCoefAndRemoveIndex(int eqnNumber, int colIndex, double& coef);

   /** Reset all coefficients to 0.0 */
   void resetCoefs();

   /** Add an equation, with associated column-indices but not coefficients.*/
   int addIndices(int eqnNumber, const int* indices, int len);

   int newCoefData_, newRHSData_;

 private:
   void deleteMemory();
   int insertNewEqn(int eqn, int insertPoint);

   int internalAddEqn(int index, const double* coefs,
                       const int* indices, int len, bool accumulate);

   std::vector<int> eqnNumbers_; //list of equation-numbers

   std::vector<fei::CSVec*> eqns_;

   std::vector<int> indices_union_; //union of all equation-indices

   int numRHSs_;     //number of right-hand-side coefficients per equation
   std::vector<std::vector<double>*> rhsCoefs_; //list of vector-pointers, each 
                                          //vector is of length numRHSs_
   bool setNumRHSsCalled_;
   bool rhsCoefsAllocated_;

   std::vector<double> dummyCoefs_;
};

std::ostream& operator<<(std::ostream& os, EqnBuffer& eq);

#endif

