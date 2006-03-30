/*
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
*/
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
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
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

// ML-headers
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"

// own headers
#include "nlnml_finelevelnoxinterface.H"

// forward declarations
class Epetra_CrsMatrix;


namespace NLNML
{


class NLNML_CoarseLevelNoxInterface : public virtual NOX::Epetra::Interface::Required
{
public:

  NLNML_CoarseLevelNoxInterface();
  void recreate();
  
  virtual ~NLNML_CoarseLevelNoxInterface();
  
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag);

  void restrict_fine_to_this();
  void prolong_this_to_fine();
  void restrict_to_next_coarser_level(); 
  void prolong_to_this_level();                                  
  
  void setP();    

  void setModifiedSystem();

  bool isFAS() { return isFASmodfied_; }

  void ApplyAllConstraints();

  int Level() { return level_; }

  void BlockMap();

private:
  NLNML_CoarseLevelNoxInterface(NLNML_CoarseLevelNoxInterface& old);
  NLNML_CoarseLevelNoxInterface operator = (const NLNML_CoarseLevelNoxInterface& old);

protected:

  int  level_;
  bool isFASmodfied_; 
  

}; // class NLNML_CoarseLevelNoxInterface
} // namespace NLNML
#endif
#endif
