// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>
#include "ml_nox_coarseinterface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface(
                                    ML_NOX::Ml_Nox_Fineinterface& interface,
                                    int level,
                                    int plevel,
                                    Epetra_CrsMatrix** P,
                                    const Epetra_BlockMap* this_RowMap,
                                    int nlevel) 
: fineinterface_(interface)
{
  if (P==0)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
         << "**ERR**: Epetra_CrsMatrix** P is 0 in constructor on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  int i;
  for (i=1; i<=level; i++)
  {
     if (P[i]==0)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
             << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }
  if (!this_RowMap)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
         << "**ERR**: Epetra_BlockMap* this_RowMap is 0 in constructor on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // store some data
  level_         = level;   // my level (ML_INCREASING only!!!)
  nlevel_        = nlevel;  // total number of levels
  ml_printlevel_ = plevel;  // set printing level
  nFcalls_       = 0;       // number of cals to the computeF function
  P_             = P;       // ptr to the vector of prolongators
  fbar_          = 0;       // the ptr to one FAS-vector
  fxbar_         = 0;       // the ptr to one FAS-vector
  isFASmodfied_  = false;
  this_RowMap_   = new Epetra_BlockMap(*this_RowMap);

  // create a series of working vectors to use in prolongations and restrictions
  if (level)
  {
    wvec_.resize(level+1);
    for (int i=0; i<level; ++i)
      wvec_[i] = new Epetra_Vector(P_[i+1]->OperatorRangeMap(),false);
    wvec_[level] = new Epetra_Vector(P_[level]->OperatorDomainMap(),false);
  }
  else
    wvec_.clear();

  return;
}

/*----------------------------------------------------------------------*
 |  recreate the coarse interface (public)                   m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::recreate(int plevel,Epetra_CrsMatrix** P, 
                                                   const Epetra_BlockMap* this_RowMap) 
{
  if (P==0)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n";
    cout << "**ERR**: Epetra_CrsMatrix** P is 0 in constructor on level " << level_ << "\n";
    cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  int i;
  for (i=1; i<=level_; i++)
  {
     if (P[i]==0)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
             << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }
  // store some data
  ml_printlevel_ = plevel;  // set printing level
  nFcalls_       = 0;       // number of cals to the computeF function
  P_             = P;       // the new prolongators
  if (!this_RowMap)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::recreate:\n"
         << "**ERR**: Epetra_BlockMap* this_RowMap is 0 on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  if (this_RowMap_) delete this_RowMap_;
  this_RowMap_   = new Epetra_BlockMap(*this_RowMap);

  // create a series of working vectors to use in prolongations and restrictions
  if (level_)
  {
    // delete the old working vectors
    for (int i=0; i<(int)wvec_.size(); ++i)
      if (wvec_[i]) delete wvec_[i];
    // create new ones
    wvec_.resize(level_+1);
    for (int i=0; i<level_; ++i)
      wvec_[i] = new Epetra_Vector(P_[i+1]->OperatorRangeMap(),false);
    wvec_[level_] = new Epetra_Vector(P_[level_]->OperatorDomainMap(),false);
  }
  else
    wvec_.clear();

  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Nox_CoarseProblem_Interface::~Nox_CoarseProblem_Interface()
{ 
   // destroy the P-hierarchy, if it still exists
   destroyP();
   if (this_RowMap_) delete this_RowMap_;
   this_RowMap_ = 0;

   for (int i=0; i<(int)wvec_.size(); ++i) 
   {
     delete wvec_[i];
     wvec_[i] = NULL;
   }
   wvec_.clear();

   return; 
}

/*----------------------------------------------------------------------*
 |  restrict a vector from level current to                  m.gee 01/05|
 |  next coarser level next                                             |
 |  returns ptr to coarse Epetra_Vector                                 |
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level(
                                                                            Epetra_Vector* thisvec, 
                                                                            int current, int next)
{
   // create a thisvec that matches the map of the restriction operator
   Epetra_Vector* xfineP = 0;
   if (current==0) // restrict from level 0 to level 1 is different from rest
   {
      // Note that the GIDs in xfine match those of the fineinterface and
      // might be different from those in P_[1]->RowMap.
      // The LIDs and the map match, so we have to copy xfine to xfineP
      // using LIDs
      xfineP = new Epetra_Vector(P_[1]->RowMap(),false);
      if (thisvec->MyLength() != xfineP->MyLength() || thisvec->GlobalLength() != xfineP->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
               << "**ERR**: mismatch in dimension of thisvec and xfineP\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      const int mylength = thisvec->MyLength();
      for (int i=0; i<mylength; i++)
         (*xfineP)[i] = (*thisvec)[i];      
   }
   else
   {
      xfineP = thisvec;
   }
   
   // Allocate matching coarse vector
   Epetra_Vector* cvec = new Epetra_Vector(P_[next]->OperatorDomainMap(),false);
   
   // restrict (transposed mutliply with prolongation)
   P_[next]->Multiply(true,*xfineP,*cvec);
   
   if (current==0)
      delete xfineP;

   return cvec;
}

/*----------------------------------------------------------------------*
 |  prolong a vector from level current to                   m.gee 01/05|
 |  next coarser level next                                             |
 |  returns ptr to coarse Epetra_Vector                                 |
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::prolong_to_this_level(Epetra_Vector* coarsevec, 
                                                                            int current, int next)
{
   // Allocate matching fine vector in ML-global indices
   Epetra_Vector* fvec = new Epetra_Vector(P_[next]->OperatorRangeMap(),false);
   
   // prolongate
   P_[next]->Multiply(false,*coarsevec,*fvec);
   
   // when prolongating to the finest level, global indexing there is different
   if (current==0)
   {
      Epetra_Vector* fine = new Epetra_Vector(fineinterface_.getMap(),false);
      if (fvec->MyLength() != fine->MyLength() || fvec->GlobalLength() != fine->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
               << "**ERR**: mismatch in dimension of fvec and fine\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      const int mylength = fvec->MyLength();
      for (int i=0; i<mylength; i++)
         (*fine)[i] = (*fvec)[i];
      delete fvec; fvec = 0;
      return fine;  
   }
   else
      return fvec;

   return NULL;
}

/*----------------------------------------------------------------------*
 |  restrict a fine level vector to this level               m.gee 12/04|
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::restrict_fine_to_this(
                                               const Epetra_Vector& xfine)
{
   int i;
   // on the finest level, just copy the vector
   if (level_==0)
   {
      Epetra_Vector* xcoarse = new Epetra_Vector(xfine);
      return xcoarse;
   }
   else // not on finest level
   {
      // Note that the GIDs in xfine match those of the fineinterface and
      // might be different from those in P_[1]->OperatorRangeMap().
      // The LIDs and the map match, so we have to copy xfine to xfineP
      // using LIDs
      Epetra_Vector* xfineP = wvec_[0]; // RangeMap of P_[1]
      if (xfine.MyLength() != xfineP->MyLength() || xfine.GlobalLength() != xfineP->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
               << "**ERR**: mismatch in dimension of xfine and xfineP\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }

      const int mylength = xfine.MyLength();
      for (i=0; i<mylength; i++)
         (*xfineP)[i] = xfine[i];      

      // loop from the finest level to this level and
      // apply series of restrictions (that is transposed prolongations)
      Epetra_Vector* fvec = xfineP;
      for (i=0; i<level_-1; i++)
      {
        Epetra_Vector* cvec = wvec_[i+1];
        P_[i+1]->Multiply(true,*fvec,*cvec);
        fvec = cvec;
      }
      Epetra_Vector* out = new Epetra_Vector(P_[level_]->OperatorDomainMap(),false); 
      P_[level_]->Multiply(true,*fvec,*out);
      return out;   
   }
   return NULL;
}

/*----------------------------------------------------------------------*
 |  prolong a this coarse level vector to the finest level   m.gee 12/04|
 |  NOTE: the calling routine is responsible of deleteing the           |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::prolong_this_to_fine(
                                            const Epetra_Vector& xcoarse)
{
   int i;
   // on the finest level, just copy the vector
   if (level_==0)
   {
      Epetra_Vector* xfine = new Epetra_Vector(xcoarse);
      return xfine;
   }
   else // we are not on the finest level
   {
      // loop from this coarsest level to the finest one
      // apply series of prolongations
      Epetra_Vector* cvec = const_cast<Epetra_Vector*>(&xcoarse);
      for (i=level_; i>0; i--)
      {
         // allocate a vector matching level i-1
         Epetra_Vector* fvec = wvec_[i-1];
         // multiply
         P_[i]->Multiply(false,*cvec,*fvec);
         cvec = fvec;
      }
      // Note that the GIDs in cvec do NOT match those of the fineinterface as
      // they match the P_[1]->RangeMap.
      // The LIDs match, so we have to copy cvec to xfine_fineinterface
      // using LIDs
      Epetra_Vector* xfine_fineinterface = new Epetra_Vector(fineinterface_.getGraph()->RowMap(),false);
      if (cvec->MyLength() != xfine_fineinterface->MyLength() ||
          cvec->GlobalLength() != xfine_fineinterface->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::prolong_this_to_fine:\n"
               << "**ERR**: mismatch in dimension of cvec and xfine_fineinterface\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      const int mylength = cvec->MyLength();
      for (i=0; i<mylength; i++)
         (*xfine_fineinterface)[i] = (*cvec)[i];
      return xfine_fineinterface;
   }
   return NULL;
}

/*----------------------------------------------------------------------*
 |  evaluate nonlinear function (public, derived)            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computeF(const Epetra_Vector& x, 
                                                   Epetra_Vector& F, 
                                                   const FillType fillFlag)
{
  bool err;
  //int  ierr;
  ++nFcalls_;

  if (level_==0)
  {
     // call fine level interface
     // in case the fillFlag is residual, sierra will not enforce constraints on it
     // As we want (need) them enforced all the time in the preconditioner
     // we switch the type to Prec
     // We leave all other types unchanged as sierra does certain things on other types 
     // which we still want to happen
     
     FillType type = fillFlag;
     //if (fillFlag == NOX::EpetraNew::Interface::Required::Residual)
       //type = NOX::EpetraNew::Interface::Required::Prec;
     err = fineinterface_.computeF(x,F,type);
     if (err==false)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: call to fine-userinterface returned false on level " << level_ << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }
  else // level_ > 0
  {
     // create Ffine and xfine matching the fine interface
     Epetra_Vector* Ffine = new Epetra_Vector(fineinterface_.getGraph()->RowMap(),false);
     Epetra_Vector* xfine = prolong_this_to_fine(x);
     
     // call fine level interface
     // in case the fillFlag is residual, sierra will not enforce constraints on it
     // As we want (need) them enforced all the time in the preconditioner
     // we switch the type to Prec
     // We leave all other types unchanged as sierra does certain things on other types 
     // which we still want to happen
     FillType type = fillFlag;
     //if (fillFlag == NOX::EpetraNew::Interface::Required::Residual)
       //type = NOX::EpetraNew::Interface::Required::Prec;
     err = fineinterface_.computeF(*xfine,*Ffine,type);
     if (err==false)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: call to fine-userinterface returned false on level " << level_ << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
     delete xfine;
     
     // get the answer Ffine back to this coarse level
     Epetra_Vector* Fcoarse = restrict_fine_to_this(*Ffine);
     delete Ffine;
     
     F.Scale(1.0,*Fcoarse);
     delete Fcoarse;
  } // level_ > 0

  // check for FAS option 
  // solving for modified system: F_c ( x_c ) - F_c ( R x_f ) + R F_f ( x_f ) -> 0 
  if (isFAS()==true)
     F.Update(-1.0,*fxbar_,1.0,*fbar_,1.0);

  return err;
}

/*----------------------------------------------------------------------*
 |  apply constraints                    (public)            m.gee 11/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Nox_CoarseProblem_Interface::ApplyAllConstraints(Epetra_Vector& gradient)
{
  if (level_==0)
  {
    // cout << "Nox_CoarseProblem_Interface::ApplyAllConstraints: fine \n"; fflush(stdout);
    fineinterface_.ApplyAllConstraints(gradient,0);
    return;
  }
  else
  {
    // cout << "Nox_CoarseProblem_Interface::ApplyAllConstraints: coarse \n"; fflush(stdout);
    Epetra_Vector* gradientfine = prolong_this_to_fine(gradient);
    fineinterface_.ApplyAllConstraints(*gradientfine,Level());
    Epetra_Vector* gradientcoarse = restrict_fine_to_this(*gradientfine);
    delete gradientfine;
    gradient.Scale(1.0,*gradientcoarse);
    delete gradientcoarse;
    return;
  }
  return;
}
/*----------------------------------------------------------------------*
 |  evaluate Jacobian           (public, derived)            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computeJacobian(const Epetra_Vector& x)
{
  cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeJacobian(...):\n"
       << "**ERR**: this  is NOT supposed to be called????????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return true;
}

/*----------------------------------------------------------------------*
 |  compute user supplied preconditioner (public, derived)   m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computePreconditioner(
                                                    const Epetra_Vector& x,
		                                    NOX::Parameter::List* precParams)
{
  cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computePreconditioner(...):\n"
       << "**ERR**: this  is NOT supposed to be called????????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return(true);
}

//-----------------------------------------------------------------------------

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
