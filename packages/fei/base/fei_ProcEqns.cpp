/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>

#include <fei_TemplateUtils.hpp>
#include <fei_ProcEqns.hpp>

//==============================================================================
ProcEqns::ProcEqns()
 : procs_(),
   eqnsPerProc_(),
   procEqnNumbers_(),
   procEqnLengths_()
{
}

//==============================================================================
ProcEqns::~ProcEqns() {
   deleteMemory();
}

//==============================================================================
ProcEqns* ProcEqns::deepCopy()
{
//
//This function returns a DEEP copy (including all data) of 'this' object.
//
   ProcEqns* dest = new ProcEqns;

   int len = procs_.size();

   dest->procs_.resize(len);
   dest->eqnsPerProc_.resize(len);

   dest->procEqnNumbers_.resize(len);
   dest->procEqnLengths_.resize(len);

   for(int i=0; i<len; i++) {
      dest->procs_[i] = procs_[i];
      dest->eqnsPerProc_[i] = eqnsPerProc_[i];

      dest->procEqnNumbers_[i] = new std::vector<int>(*(procEqnNumbers_[i]));
      dest->procEqnLengths_[i] = new std::vector<int>(*(procEqnLengths_[i]));
   }

   return(dest);
}

//==============================================================================
void ProcEqns::deleteMemory() {
  for(unsigned i=0; i<procEqnNumbers_.size(); i++) {
    delete procEqnNumbers_[i];
    delete procEqnLengths_[i];
  }
}

//==============================================================================
void ProcEqns::addEqn(int eqnNumber, int proc) {

   internalAddEqn(eqnNumber, 0, proc);
}

//==============================================================================
void ProcEqns::addEqn(int eqnNumber, int eqnLength, int proc) {
   internalAddEqn(eqnNumber, eqnLength, proc);
}

//==============================================================================
void ProcEqns::internalAddEqn(int eqnNumber, int eqnLength, int proc) {
//
//This function adds proc to the recvProcs_ list if it isn't already
//present, and adds eqnNumber to the correct row of the procEqnNumbers_
//table if eqnNumber isn't already in there.
//

  //is proc already in our list of procs?
  std::vector<int>::iterator
    p_iter = std::lower_bound(procs_.begin(),procs_.end(), proc);

  unsigned offset = p_iter - procs_.begin();

  if (p_iter == procs_.end() || proc != *p_iter) {
    //proc was NOT already present, so
    //we need to insert it, and insert new rows in the tables
    //procEqnNumbers_ and procEqnLengths.

    procs_.insert(p_iter, proc);

    procEqnNumbers_.insert(procEqnNumbers_.begin()+offset,new std::vector<int>(1));
    (*(procEqnNumbers_[offset]))[0] = eqnNumber;

    procEqnLengths_.insert(procEqnLengths_.begin()+offset, new std::vector<int>(1));
    (*(procEqnLengths_[offset]))[0] = eqnLength;

    eqnsPerProc_.insert(eqnsPerProc_.begin()+offset, 1);
  }
  else {
    //proc was already in the procs_ list.
    //is eqnNumber already in our list of eqns for proc?
    //if not, add it.

    std::vector<int>& procEqnNums = *(procEqnNumbers_[offset]);
    std::vector<int>::iterator pe_iter =
      std::lower_bound(procEqnNums.begin(),procEqnNums.end(), eqnNumber);

    unsigned offset2 = pe_iter - procEqnNums.begin();

    if (pe_iter == procEqnNums.end() || eqnNumber != *pe_iter) {
      procEqnNumbers_[offset]->insert(procEqnNumbers_[offset]->begin()+offset2,eqnNumber);
      procEqnLengths_[offset]->insert(procEqnLengths_[offset]->begin()+offset2,eqnLength);
      eqnsPerProc_[offset] = procEqnNumbers_[offset]->size();
    }
    else {
      (*(procEqnLengths_[offset]))[offset2] = eqnLength;
    }
  }
}

//==============================================================================
void ProcEqns::setProcEqnLengths(int* eqnNumbers, int* eqnLengths, int len)
{
   if (len == 0) return;

   int numProcs = procs_.size();

   for(int i=0; i<numProcs; i++) {
      for(int j=0; j<eqnsPerProc_[i]; j++) {
         int eqn_j = (*(procEqnNumbers_[i]))[j];
         int ins = -1;
	 int index = fei::binarySearch( eqn_j,
					    eqnNumbers, len, ins);

         if (index < 0) {
            fei::console_out() << "fei: ProcEqns::setProcEqnLengths: ERROR, "
                 << eqn_j << " not found." << FEI_ENDL;
            std::abort();
         }

         (*(procEqnLengths_[i]))[j] = eqnLengths[index];
      }
   }
}

