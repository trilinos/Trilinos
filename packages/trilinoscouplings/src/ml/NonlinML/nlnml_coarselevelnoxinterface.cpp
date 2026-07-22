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
#include "nlnml_coarselevelnoxinterface.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_CoarseLevelNoxInterface::NLNML_CoarseLevelNoxInterface(
                 RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface,
                 int level,int outlevel,
                 RefCountPtr< vector< RefCountPtr<Epetra_CrsMatrix> > > P,
                 const Epetra_BlockMap& this_bmap,
                 int maxlevel) :
level_(level),
maxlevel_(maxlevel),
outputlevel_(outlevel),
isFASmodfied_(false),
nFcalls_(0),
fineinterface_(finterface),
P_(P),
fbar_(NULL),
fxbar_(NULL)
{
  this_bmap_ = rcp(new Epetra_BlockMap(this_bmap));

  // create a series of working vectors to be used on prolongation/restriction
  wvec_.clear();
  if (Level())
  {
    wvec_.resize(Level()+1);
    for (int i=0; i<Level(); ++i)
      wvec_[i] = rcp(new Epetra_Vector((*P)[i+1]->OperatorRangeMap(),false));
    wvec_[Level()] = rcp(new Epetra_Vector((*P)[Level()]->OperatorDomainMap(),false));
  }



  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::recreate()
{

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_CoarseLevelNoxInterface::~NLNML_CoarseLevelNoxInterface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nonlinear function (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_CoarseLevelNoxInterface::computeF(
                                 const Epetra_Vector& x, Epetra_Vector& F,
               const FillType fillFlag)
{
  bool err;
  if (!Level())
  {
    err = fineinterface_->computeF(x,F,fillFlag);
    if (err==false)
    {
      cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::computeF:\n"
           << "**ERR**: call to fine-userinterface returned false on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
  }
  else
  {
    RefCountPtr<Epetra_Vector> Ffine = rcp(new Epetra_Vector(fineinterface_->getGraph()->RowMap(),false));
    RefCountPtr<Epetra_Vector> xfine = rcp(prolong_this_to_fine(x));
    err = fineinterface_->computeF(*xfine,*Ffine,fillFlag);
    if (err==false)
    {
      cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::computeF:\n"
           << "**ERR**: call to fine-userinterface returned false on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    RefCountPtr<Epetra_Vector> Fcoarse = rcp(restrict_fine_to_this(*Ffine));
    F.Scale(1.0,*Fcoarse);
  }
  if (isFAS())
    F.Update(-1.0,*fxbar_,1.0,*fbar_,1.0);
  return err;
}

/*----------------------------------------------------------------------*
 |  restrict from fine to this level (public)                 m.gee 3/06|
 *----------------------------------------------------------------------*/
Epetra_Vector*  NLNML::NLNML_CoarseLevelNoxInterface::restrict_fine_to_this(
                                                 const Epetra_Vector& xfine)
{
  if (!Level())
  {
    Epetra_Vector* xcoarse = new Epetra_Vector(xfine);
    return xcoarse;
  }
  else
  {
    // Note that the GIDs in xfine match those of the fineinterface and
    // might be different from those in P_[1]->OperatorRangeMap().
    // The LIDs and the map match, so we have to copy xfine to xfineP
    // using LIDs
    Epetra_Vector* xfineP = wvec_[0].get(); // RangeMap of P_[1]
    if (xfine.MyLength() != xfineP->MyLength() || xfine.GlobalLength() != xfineP->GlobalLength())
    {
        cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::restrict_fine_to_this:\n"
             << "**ERR**: mismatch in dimension of xfine and xfineP\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    const int mylength = xfine.MyLength();
    for (int i=0; i<mylength; i++)
      (*xfineP)[i] = xfine[i];

    // loop from the finest level to this level and
    // apply series of restrictions (that is transposed prolongations)
    Epetra_Vector* fvec = xfineP;
    for (int i=0; i<Level()-1; i++)
    {
      Epetra_Vector* cvec = wvec_[i+1].get();
      (*P_)[i+1]->Multiply(true,*fvec,*cvec);
      fvec = cvec;
    }
    Epetra_Vector* out = new Epetra_Vector((*P_)[Level()]->OperatorDomainMap(),false);
    (*P_)[Level()]->Multiply(true,*fvec,*out);
    return out;
  }
}

/*----------------------------------------------------------------------*
 |  prolongate from this level to fine (public)               m.gee 3/06|
 *----------------------------------------------------------------------*/
Epetra_Vector* NLNML::NLNML_CoarseLevelNoxInterface::prolong_this_to_fine(
                                                const Epetra_Vector& xcoarse)
{
  if (!Level())
    return new Epetra_Vector(xcoarse);
  else
  {
    Epetra_Vector* cvec = const_cast<Epetra_Vector*>(&xcoarse);
    for (int i=Level(); i>0; i--)
    {
       // allocate a vector matching level i-1
       Epetra_Vector* fvec = wvec_[i-1].get();
       // multiply
       (*P_)[i]->Multiply(false,*cvec,*fvec);
       cvec = fvec;
    }
    // Note that the GIDs in cvec do NOT match those of the fineinterface as
    // they match the P_[1]->RangeMap.
    // The LIDs match, so we have to copy cvec to xfine_fineinterface
    // using LIDs
    Epetra_Vector* xfine_fineinterface = new Epetra_Vector(fineinterface_->getGraph()->RowMap(),false);
    if (cvec->MyLength() != xfine_fineinterface->MyLength() ||
        cvec->GlobalLength() != xfine_fineinterface->GlobalLength())
    {
        cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::prolong_this_to_fine:\n"
             << "**ERR**: mismatch in dimension of cvec and xfine_fineinterface\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    const int mylength = cvec->MyLength();
    for (int i=0; i<mylength; i++)
       (*xfine_fineinterface)[i] = (*cvec)[i];
    return xfine_fineinterface;
  }
}

/*----------------------------------------------------------------------*
 |  restrict from this to next coarser level (public)         m.gee 3/06|
 *----------------------------------------------------------------------*/
Epetra_Vector* NLNML::NLNML_CoarseLevelNoxInterface::restrict_to_next_coarser_level(
                                                  Epetra_Vector* thisvec,
                                                  int current, int next)
{
  Epetra_Vector* xfineP = 0;
  if (current==0)
  {
    xfineP = new Epetra_Vector((*P_)[1]->RowMap(),false);
    if (thisvec->MyLength() != xfineP->MyLength() || thisvec->GlobalLength() != xfineP->GlobalLength())
    {
      cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::restrict_to_next_coarser_level:\n"
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
  Epetra_Vector* cvec = new Epetra_Vector((*P_)[next]->OperatorDomainMap(),false);
  (*P_)[next]->Multiply(true,*xfineP,*cvec);
  if (current==0) delete xfineP;
  return cvec;
}


/*----------------------------------------------------------------------*
 |  prolongate from next coarser level to this level (public) m.gee 3/06|
 *----------------------------------------------------------------------*/
Epetra_Vector* NLNML::NLNML_CoarseLevelNoxInterface::prolong_to_this_level(
                                                  Epetra_Vector* coarsevec,
                                                  int current, int next)
{
  Epetra_Vector* fvec = new Epetra_Vector((*P_)[next]->OperatorRangeMap(),false);
  (*P_)[next]->Multiply(false,*coarsevec,*fvec);
  if (current==0)
  {
    Epetra_Vector* fine = new Epetra_Vector(fineinterface_->getMap(),false);
    if (fvec->MyLength() != fine->MyLength() || fvec->GlobalLength() != fine->GlobalLength())
    {
      cout << "**ERR**: NLNML::NLNML_CoarseLevelNoxInterface::prolong_to_this_level:\n"
           << "**ERR**: mismatch in dimension of fvec and fine\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    const int mylength = fvec->MyLength();
    for (int i=0; i<mylength; ++i)
      (*fine)[i] = (*fvec)[i];
    delete fvec;
    return fine;
  }
  else
    return fvec;
}

/*----------------------------------------------------------------------*
 |  set ptr to all prolongation operators (public)            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::setP()
{
  return;
}


/*----------------------------------------------------------------------*
 |  make application apply all constraints (public)           m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::ApplyAllConstraints(
                                                  Epetra_Vector& gradient)
{
  if (!Level())
    fineinterface_->ApplyAllConstraints(gradient,0);
  else
  {
    RefCountPtr<Epetra_Vector> gradientfine = rcp(prolong_this_to_fine(gradient));
    fineinterface_->ApplyAllConstraints(*gradientfine,Level());
    RefCountPtr<Epetra_Vector> gradientcoarse = rcp(restrict_fine_to_this(*gradientfine));
    gradient.Scale(1.0,*gradientcoarse);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  return this levels blockmap (public)                      m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::BlockMap()
{
  return;
}









#endif
