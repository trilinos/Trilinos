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
/*!
 * \file nlnml_finelevelnoxinterface.H
 *
 * \class NLNML_finelevelnoxinterface
 *
 * \brief class NLNML_coarselevelnoxinterface used by the nonlinear preconditioner
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "ml_common.h"
#include "ml_include.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
#ifndef NLNML_COARSELEVELNOXINTERFACE_H
#define NLNML_COARSELEVELNOXINTERFACE_H

#include <iostream>

// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX headers
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

#include "Teuchos_RefCountPtr.hpp"

// ML-headers
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"

// own headers
#include "nlnml_finelevelnoxinterface.H"

// forward declarations
class Epetra_CrsMatrix;

using namespace Teuchos;
namespace NLNML
{


class NLNML_CoarseLevelNoxInterface : public virtual NOX::Epetra::Interface::Required
{
public:

  NLNML_CoarseLevelNoxInterface(
                     RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface,
                     int level,int outlevel,
                     RefCountPtr< std::vector< RefCountPtr<Epetra_CrsMatrix> > > P,
                     const Epetra_BlockMap& this_bmap,
                     int maxlevel);
  void recreate();

  virtual ~NLNML_CoarseLevelNoxInterface();

  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
      const FillType fillFlag);

  Epetra_Vector* restrict_fine_to_this(const Epetra_Vector& xfine);
  Epetra_Vector* prolong_this_to_fine(const Epetra_Vector& xcoarse);
  Epetra_Vector* restrict_to_next_coarser_level(Epetra_Vector* thisvec, int current, int next);
  Epetra_Vector* prolong_to_this_level(Epetra_Vector* coarsevec, int current, int next);

  void setP();

  bool setModifiedSystem(bool flag, Epetra_Vector* fbar, Epetra_Vector* fxbar)
  { isFASmodfied_ = flag; fbar_ = fbar; fxbar_ = fxbar; return true; }

  bool isFAS() { return isFASmodfied_; }

  void ApplyAllConstraints(Epetra_Vector& gradient);

  int Level() { return level_; }

  void BlockMap();

  bool resetnumcallscomputeF() { nFcalls_=0; return true;}

private:
  NLNML_CoarseLevelNoxInterface(NLNML_CoarseLevelNoxInterface& old);
  NLNML_CoarseLevelNoxInterface operator = (const NLNML_CoarseLevelNoxInterface& old);

  inline int Output() const { return outputlevel_; }
  inline int Level() const { return level_; }

protected:

  int                                                    level_;
  int                                                    maxlevel_;
  int                                                    outputlevel_;
  bool                                                   isFASmodfied_;
  int                                                    nFcalls_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface>        fineinterface_;
  RefCountPtr< std::vector< RefCountPtr<Epetra_CrsMatrix> > > P_;
  RefCountPtr<Epetra_BlockMap>                           this_bmap_;
  std::vector<RefCountPtr<Epetra_Vector> >                    wvec_;
  Epetra_Vector*                                         fbar_;
  Epetra_Vector*                                         fxbar_;
}; // class NLNML_CoarseLevelNoxInterface
} // namespace NLNML
#endif
#endif
