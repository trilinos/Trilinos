/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_EqnBuffer.hpp>
#include <fei_SSVec.hpp>

#include <feiArray.hpp>
#include <fei_TemplateUtils.hpp>

//==============================================================================
EqnBuffer::EqnBuffer()
 : newCoefData_(0),
   newRHSData_(0),
   eqnNumbers_(0, 8000),
   eqns_(0,8000),
   indices_union_(0, 2000),
   numRHSs_(1),
   rhsCoefs_(0,8000),
   setNumRHSsCalled_(false),
   rhsCoefsAllocated_(false),
   dummyCoefs_(0, 32),
   dummyCoefsLen_(0)
{
}

//==============================================================================
EqnBuffer::EqnBuffer(const EqnBuffer& src)
 : newCoefData_(0),
   newRHSData_(0),
   eqnNumbers_(0, 8000),
   eqns_(0,8000),
   indices_union_(0, 2000),
   numRHSs_(1),
   rhsCoefs_(0,8000),
   setNumRHSsCalled_(false),
   rhsCoefsAllocated_(false),
   dummyCoefs_(0, 32),
   dummyCoefsLen_(0)
{
  *this = src;
}

//==============================================================================
EqnBuffer& EqnBuffer::operator=(const EqnBuffer& src)
{
   int i, len = src.eqnNumbers_.length();

   eqnNumbers_ = src.eqnNumbers_;
   eqns_.resize(src.eqns_.length());

   numRHSs_ = src.numRHSs_;

   for(i=0; i<len; i++) {
     //eqns_ is a table. Each row of the table needs to be allocated and
     //copied here. We'll use the SSVec copy constructor to copy the
     //contents of each existing row into the 'dest' rows.

      //first get a pointer to the row,
      SSVec* row = src.eqns_[i];

      //now allocate the eqns_ row and the coefs row
      eqns_[i] = new SSVec(*row);

      //since we allow for multiple rhs's, rhsCoefs_ is a table too...
      feiArray<double>* rhsCoefs = src.rhsCoefs_[i];

      rhsCoefs_.append( new feiArray<double>(*rhsCoefs) );
   }

   return(*this);
}

//==============================================================================
EqnBuffer::~EqnBuffer() {
   deleteMemory();
}

//==============================================================================
EqnBuffer* EqnBuffer::deepCopy()
{
   EqnBuffer* dest = new EqnBuffer;

   *dest = *this;
 
   return(dest);
}

//==============================================================================
void EqnBuffer::deleteMemory() {
   for(int i=0; i<getNumEqns(); i++) {
      delete eqns_[i];

      delete rhsCoefs_[i];
   }

   eqns_.reAllocate(0);
   rhsCoefs_.reAllocate(0);
   numRHSs_ = 0;
}

//==============================================================================
int EqnBuffer::getEqnIndex(int eqn) {

   return(snl_fei::binarySearch(eqn, eqnNumbers_));
}

//==============================================================================
void EqnBuffer::setNumRHSs(int n) {

   if (n <= 0) { return;}

   numRHSs_ = n;

   if (getNumEqns() <= 0) return;

   for(int i=0; i<getNumEqns(); i++) {
     feiArray<double>* rhsCoefs = rhsCoefs_[i];
     rhsCoefs->resize(numRHSs_);

     (*rhsCoefs) = 0.0; // uses feiArray operator= for initialization
   }
}

//==============================================================================
int EqnBuffer::addRHS(int eqnNumber, int rhsIndex, double value,
		      bool accumulate)
{
   int index = snl_fei::binarySearch(eqnNumber, eqnNumbers_);

   if (index < 0) {
      FEI_CERR << "(deep in FEI) EqnBuffer::addRHS: ERROR, eqnNumber " << eqnNumber
           << " not found in send eqns." << FEI_ENDL;
      return(-1);
   }

   feiArray<double>* rhsCoefs = rhsCoefs_[index];

   if ( rhsCoefs->length() <= rhsIndex) setNumRHSs(rhsIndex+1);

   if (accumulate==true) (*rhsCoefs)[rhsIndex] += value;
   else (*rhsCoefs)[rhsIndex] = value;

   return(0);
}
//==============================================================================
int EqnBuffer::isInIndices(int eqn)
{
  //
  //This function checks the indices_ table to see if 'eqn' is present.
  //If it is, the appropriate row index into the table is returned.
  //-1 is return otherwise.
  //
  if (indices_union_.length() > 0) {
    int index = snl_fei::binarySearch(eqn, indices_union_.dataPtr(),
                                      indices_union_.length());
    if (index < 0) return(-1);
  }

  int numEqns = getNumEqns(), index;
  SSVec** eqnsPtr = eqns_.dataPtr();
  for(int i=0; i<numEqns; i++) {
    feiArray<int>& indices = eqnsPtr[i]->indices();
    index = snl_fei::binarySearch(eqn, indices.dataPtr(),
                                  indices.length());
    if (index > -1) return(i);
  }

  return(-1);
}

//==============================================================================
int EqnBuffer::addEqn(int eqnNumber, const double* coefs, const int* indices,
                       int len, bool accumulate, bool create_indices_union) 
{
  if (len <= 0) return(0);

  int err, insertPoint = -1;
  int index = snl_fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

  if (index < 0) {
    //if eqnNumber isn't already present, insert a new entry into the
    //appropriate data structures.
    err = insertNewEqn(eqnNumber, insertPoint);
    if (err) {return(err);}
    index = insertPoint;
  }

  //Now add the coef/index values.
  err = internalAddEqn(index, coefs, indices, len, accumulate);

  if (create_indices_union) {
    for(int i=0; i<len; ++i) {
      snl_fei::sortedListInsert(indices[i], indices_union_);
    }
  }

  return(err);
}

//==============================================================================
int EqnBuffer::getCoef(int eqnNumber, int colIndex, double& coef)
{
  int eqnLoc = snl_fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = snl_fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(-1);

  coef = eqns_[eqnLoc]->coefs()[colLoc];
  return(0);
}

//==============================================================================
int EqnBuffer::removeIndex(int eqnNumber, int colIndex)
{
  int eqnLoc = snl_fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = snl_fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(0);

  feiArray<int>& indices = eqns_[eqnLoc]->indices();
  feiArray<double>& coefs= eqns_[eqnLoc]->coefs();

  int len = indices.length();

  int* indPtr = indices.dataPtr();
  double* coefPtr = coefs.dataPtr();

  for(int i=len-1; i>colLoc; --i) {
    indPtr[i-1] = indPtr[i];
    coefPtr[i-1] = coefPtr[i];
  }

  indices.resize(len-1);
  coefs.resize(len-1);

  return(0);
}

//==============================================================================
int EqnBuffer::getCoefAndRemoveIndex(int eqnNumber, int colIndex, double& coef)
{
  int eqnLoc = snl_fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = snl_fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(-1);

  feiArray<int>& indices = eqns_[eqnLoc]->indices();
  feiArray<double>& coefs= eqns_[eqnLoc]->coefs();

  coef = coefs[colLoc];
  int len = indices.length();

  int* indPtr = indices.dataPtr();
  double* coefPtr = coefs.dataPtr();

  for(int i=len-1; i>colLoc; --i) {
    indPtr[i-1] = indPtr[i];
    coefPtr[i-1] = coefPtr[i];
  }

  indices.resize(len-1);
  coefs.resize(len-1);

  return(0);
}

//==============================================================================
int EqnBuffer::addEqns(EqnBuffer& inputEqns, bool accumulate)
{
  int* eqnNums = inputEqns.eqnNumbersPtr().dataPtr();
  SSVec** eqs = inputEqns.eqns().dataPtr();

  int numRHSs = inputEqns.getNumRHSs();
  feiArray<double>** rhsCoefs = inputEqns.rhsCoefsPtr()->dataPtr();

  for(int i=0; i<inputEqns.getNumEqns(); i++) {
    feiArray<int>& indices_i  = eqs[i]->indices();
    feiArray<double>& coefs_i = eqs[i]->coefs();

    int err = addEqn(eqnNums[i], coefs_i.dataPtr(), indices_i.dataPtr(),
		    eqs[i]->length(), accumulate);
    if (err) return(err);

    if (numRHSs > 0) {
      for(int j=0; j<numRHSs; ++j) {
	addRHS(eqnNums[i], j, (*(rhsCoefs[i]))[j], accumulate);
      }
    }
  }

  return(0);
}

//==============================================================================
int EqnBuffer::insertNewEqn(int eqn, int insertPoint)
{
  //private function. We can safely assume that insertPoint is the correct
  //offset at which to insert the new equation.
  try {
    eqnNumbers_.insert(eqn, insertPoint);

    SSVec* newEqn = new SSVec;
    eqns_.insert(newEqn, insertPoint);

    if (numRHSs_ <= 0) return(-1);

    feiArray<double>* newRhsCoefRow = new feiArray<double>(numRHSs_, 1);
    *newRhsCoefRow = 0.0;
    rhsCoefs_.insert(newRhsCoefRow, insertPoint);
  }
  catch (std::runtime_error& exc) {
    FEI_CERR << exc.what() << FEI_ENDL;
    return(-1);
  }

  return(0);
}

//==============================================================================
int EqnBuffer::internalAddEqn(int index, const double* coefs,
			      const int* indices, int len, bool accumulate) 
{
  //
  //Private EqnBuffer function. We can safely assume that this function is only
  //called if indices_ and coefs_ already contain an 'index'th row.
  //

  SSVec& eqn = *(eqns_[index]);
  int err = 0;

  if (accumulate) {
    err = eqn.addEntries(len, coefs, indices);
  }
  else {
    err = eqn.putEntries(len, coefs, indices);
  }

  return(err);
}

//==============================================================================
void EqnBuffer::resetCoefs() {

   for(int i=0; i<getNumEqns(); i++) {
     feiArray<double>& coefRow = eqns_[i]->coefs();
     feiArray<double>* rhsCoefRow = rhsCoefs_[i];

     //now call the initialization 'operator=' for these arrays...
     coefRow = 0.0;
     *rhsCoefRow = 0.0;
   }
}

//==============================================================================
int EqnBuffer::addIndices(int eqnNumber, const int* indices, int len)
{
   int err, insertPoint = -1;
   int index = snl_fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

   //(we're adding dummy coefs as well, even though there are no
   //incoming coefs at this point).

   if (dummyCoefsLen_ < len) {
     dummyCoefs_.resize(len);
     dummyCoefs_ = 0.0;
     dummyCoefsLen_ = len;
   }

   if (index < 0) {
     //if eqnNumber was not already present, insert new equation

     err = insertNewEqn(eqnNumber, insertPoint);
     if (err) {return(err);}
     index = insertPoint;
   }

   err = internalAddEqn(index, dummyCoefs_.dataPtr(), indices, len, true);
   return(err);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, EqnBuffer& eq)
{
  feiArray<int>& eqnNums = eq.eqnNumbersPtr();
  feiArray<feiArray<double>*>& rhsCoefs = *(eq.rhsCoefsPtr());

  os << "#ereb num-eqns: " << eqnNums.length() << FEI_ENDL;
  for(int i=0; i<eqnNums.length(); i++) {
    os << "#ereb eqn " << eqnNums[i] << ": ";

    feiArray<int>& inds = eq.eqns()[i]->indices();
    feiArray<double>& cfs = eq.eqns()[i]->coefs();

    for(int j=0; j<inds.length(); j++) {
      os << "("<<inds[j] << "," << cfs[j] << ") ";
    }

    os << " rhs: ";
    feiArray<double>& rhs = *(rhsCoefs[i]);
    for(int k=0; k<rhs.length(); k++) {
      os << rhs[k] << ", ";
    }

    os << FEI_ENDL;
  }

  return(os);
}
