/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_EqnBuffer.hpp>
#include <fei_CSVec.hpp>

#include <fei_TemplateUtils.hpp>

//==============================================================================
EqnBuffer::EqnBuffer()
 : newCoefData_(0),
   newRHSData_(0),
   eqnNumbers_(0),
   eqns_(),
   indices_union_(0),
   numRHSs_(1),
   rhsCoefs_(),
   setNumRHSsCalled_(false),
   rhsCoefsAllocated_(false),
   dummyCoefs_()
{
}

//==============================================================================
EqnBuffer::EqnBuffer(const EqnBuffer& src)
 : newCoefData_(0),
   newRHSData_(0),
   eqnNumbers_(0),
   eqns_(),
   indices_union_(0),
   numRHSs_(1),
   rhsCoefs_(),
   setNumRHSsCalled_(false),
   rhsCoefsAllocated_(false),
   dummyCoefs_()
{
  *this = src;
}

//==============================================================================
EqnBuffer& EqnBuffer::operator=(const EqnBuffer& src)
{
   int i, len = src.eqnNumbers_.size();

   eqnNumbers_ = src.eqnNumbers_;
   eqns_.resize(src.eqns_.size());

   numRHSs_ = src.numRHSs_;

   for(i=0; i<len; i++) {
     //eqns_ is a table. Each row of the table needs to be allocated and
     //copied here. We'll use the fei::CSVec copy constructor to copy the
     //contents of each existing row into the 'dest' rows.

      //first get a pointer to the row,
      fei::CSVec* row = src.eqns_[i];

      //now allocate the eqns_ row and the coefs row
      eqns_[i] = new fei::CSVec(*row);

      //since we allow for multiple rhs's, rhsCoefs_ is a table too...
      std::vector<double>* rhsCoefs = src.rhsCoefs_[i];

      rhsCoefs_.push_back( new std::vector<double>(*rhsCoefs) );
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

   eqns_.clear();
   rhsCoefs_.clear();
   numRHSs_ = 0;
}

//==============================================================================
int EqnBuffer::getEqnIndex(int eqn) {

   return(fei::binarySearch(eqn, eqnNumbers_));
}

//==============================================================================
void EqnBuffer::setNumRHSs(int n) {

   if (n <= 0) { return;}

   numRHSs_ = n;

   if (getNumEqns() <= 0) return;

   for(int i=0; i<getNumEqns(); i++) {
     std::vector<double>* rhsCoefs = rhsCoefs_[i];
     rhsCoefs->assign(numRHSs_, 0.0);
   }
}

//==============================================================================
int EqnBuffer::addRHS(int eqnNumber, int rhsIndex, double value,
		      bool accumulate)
{
   int index = fei::binarySearch(eqnNumber, eqnNumbers_);

   if (index < 0) {
      fei::console_out() << "(deep in FEI) EqnBuffer::addRHS: ERROR, eqnNumber " << eqnNumber
           << " not found in send eqns." << FEI_ENDL;
      return(-1);
   }

   std::vector<double>* rhsCoefs = rhsCoefs_[index];

   if ( (int)rhsCoefs->size() <= rhsIndex) setNumRHSs(rhsIndex+1);

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
  if (indices_union_.size() > 0) {
    int index = fei::binarySearch(eqn, &indices_union_[0], indices_union_.size());
    if (index < 0) return(-1);
  }

  int numEqns = getNumEqns(), index;
  fei::CSVec** eqnsPtr = &eqns_[0];
  for(int i=0; i<numEqns; i++) {
    std::vector<int>& indices = eqnsPtr[i]->indices();
    index = fei::binarySearch(eqn, &indices[0], indices.size());
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
  int index = fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

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
      fei::sortedListInsert(indices[i], indices_union_);
    }
  }

  return(err);
}

//==============================================================================
int EqnBuffer::getCoef(int eqnNumber, int colIndex, double& coef)
{
  int eqnLoc = fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(-1);

  coef = eqns_[eqnLoc]->coefs()[colLoc];
  return(0);
}

//==============================================================================
int EqnBuffer::removeIndex(int eqnNumber, int colIndex)
{
  int eqnLoc = fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(0);

  std::vector<int>& indices = eqns_[eqnLoc]->indices();
  std::vector<double>& coefs= eqns_[eqnLoc]->coefs();

  int len = indices.size();

  int* indPtr = &indices[0];
  double* coefPtr = &coefs[0];

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
  int eqnLoc = fei::binarySearch(eqnNumber, eqnNumbers_);
  if (eqnLoc < 0) return(-1);

  int colLoc = fei::binarySearch(colIndex, eqns_[eqnLoc]->indices());
  if (colLoc < 0) return(-1);

  std::vector<int>& indices = eqns_[eqnLoc]->indices();
  std::vector<double>& coefs= eqns_[eqnLoc]->coefs();

  coef = coefs[colLoc];
  int len = indices.size();

  int* indPtr = &indices[0];
  double* coefPtr = &coefs[0];

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
  if (inputEqns.eqnNumbers().size() < 1) {
    return(0);
  }

  int* eqnNums = &(inputEqns.eqnNumbers()[0]);
  fei::CSVec** eqs = &(inputEqns.eqns()[0]);

  int numRHSs = inputEqns.getNumRHSs();
  std::vector<double>** rhsCoefs = &((*(inputEqns.rhsCoefsPtr()))[0]);

  for(int i=0; i<inputEqns.getNumEqns(); i++) {
    std::vector<int>& indices_i  = eqs[i]->indices();
    std::vector<double>& coefs_i = eqs[i]->coefs();

    int err = addEqn(eqnNums[i], &coefs_i[0], &indices_i[0],
		    eqs[i]->size(), accumulate);
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
    eqnNumbers_.insert(eqnNumbers_.begin()+insertPoint, eqn);

    fei::CSVec* newEqn = new fei::CSVec;
    eqns_.insert(eqns_.begin()+insertPoint, newEqn);

    if (numRHSs_ <= 0) return(-1);

    std::vector<double>* newRhsCoefRow = new std::vector<double>(numRHSs_, 0.0);
    rhsCoefs_.insert(rhsCoefs_.begin()+insertPoint, newRhsCoefRow);
  }
  catch (std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
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

  fei::CSVec& eqn = *(eqns_[index]);

  if (accumulate) {
    for(int i=0; i<len; ++i) {
      fei::add_entry(eqn, indices[i], coefs[i]);
    }
  }
  else {
    for(int i=0; i<len; ++i) {
      fei::put_entry(eqn, indices[i], coefs[i]);
    }
  }

  return(0);
}

//==============================================================================
void EqnBuffer::resetCoefs() {

  for(int i=0; i<getNumEqns(); i++) {
    fei::set_values(*eqns_[i], 0.0);
  }
}

//==============================================================================
int EqnBuffer::addIndices(int eqnNumber, const int* indices, int len)
{
   int err = 0, insertPoint = -1;
   int index = fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

   //(we're adding dummy coefs as well, even though there are no
   //incoming coefs at this point).

   if ((int)dummyCoefs_.size() < len) {
     dummyCoefs_.assign(len, 0.0);
   }

   if (index < 0) {
     //if eqnNumber was not already present, insert new equation

     err = insertNewEqn(eqnNumber, insertPoint);
     if (err) {return(err);}
     index = insertPoint;
   }

   if (len > 0) {
     err = internalAddEqn(index, &dummyCoefs_[0], indices, len, true);
   }
   return(err);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, EqnBuffer& eq)
{
  std::vector<int>& eqnNums = eq.eqnNumbers();
  std::vector<std::vector<double>*>& rhsCoefs = *(eq.rhsCoefsPtr());

  os << "#ereb num-eqns: " << eqnNums.size() << FEI_ENDL;
  for(size_t i=0; i<eqnNums.size(); i++) {
    os << "#ereb eqn " << eqnNums[i] << ": ";

    std::vector<int>& inds = eq.eqns()[i]->indices();
    std::vector<double>& cfs = eq.eqns()[i]->coefs();

    for(size_t j=0; j<inds.size(); j++) {
      os << "("<<inds[j] << "," << cfs[j] << ") ";
    }

    os << " rhs: ";
    std::vector<double>& rhs = *(rhsCoefs[i]);
    for(size_t k=0; k<rhs.size(); k++) {
      os << rhs[k] << ", ";
    }

    os << FEI_ENDL;
  }

  return(os);
}
