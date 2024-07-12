// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_FiniteDifferenceColoringWithUpdate.H"
#include "Epetra_MapColoring.h"
#include "Epetra_CrsMatrix.h"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

using namespace NOX;
using namespace NOX::Epetra;


// Constructor with no frills
FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(
      Teuchos::ParameterList& printingParams_,
      const Teuchos::RCP<Interface::Required>& i_,
      const NOX::Epetra::Vector& initialGuess_,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
      double beta_, double alpha_):
  FiniteDifference(printingParams_, i_, initialGuess_, beta_, alpha_),
  jacobianComputed(false),
  use_update(false),
  use_probing_diags(false),
  colorMap_(colorMap),
  updateColorMap_(Teuchos::null)
{
}

  //! Constructor with graph
 FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(
      Teuchos::ParameterList& printingParams_,
      const Teuchos::RCP<Interface::Required>& i_,
      const NOX::Epetra::Vector& initialGuess_,
      const Teuchos::RCP<Epetra_CrsGraph>& rawGraph_,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
      double beta_, double alpha_):
  FiniteDifference(printingParams_, i_, initialGuess_, rawGraph_, beta_, alpha_),
  jacobianComputed(false),
  use_update(false),
  use_probing_diags(false),
  colorMap_(colorMap),
  updateColorMap_(Teuchos::null)
{
}


FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(
      Teuchos::ParameterList& printingParams_,
      const Teuchos::RCP<Interface::Required>& i_,
      const NOX::Epetra::Vector& initialGuess_,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
      const Teuchos::RCP<Epetra_MapColoring>& updatecolorMap,
      double beta_, double alpha_):
  FiniteDifference(printingParams_, i_, initialGuess_, beta_, alpha_),
  jacobianComputed(false),
  use_update(true),
  use_probing_diags(false),
  colorMap_(colorMap),
  updateColorMap_(updatecolorMap)
{
}

FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(
      Teuchos::ParameterList& printingParams_,
      const Teuchos::RCP<Interface::Required>& i_,
      const NOX::Epetra::Vector& initialGuess_,
      const Teuchos::RCP<Epetra_CrsGraph>& rawGraph_,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
      const Teuchos::RCP<Epetra_MapColoring>& updatecolorMap,
      double beta_, double alpha_):
  FiniteDifference(printingParams_, i_, initialGuess_, rawGraph_, beta_, alpha_),
  jacobianComputed(false),
  use_update(true),
  use_probing_diags(false),
  colorMap_(colorMap),
  updateColorMap_(updatecolorMap)
{
}

FiniteDifferenceColoringWithUpdate::~FiniteDifferenceColoringWithUpdate(){}

/*****************************************************/
bool FiniteDifferenceColoringWithUpdate::computeJacobian(const Epetra_Vector& x){
  return computeJacobian(x, *this);
}


/*****************************************************/
bool FiniteDifferenceColoringWithUpdate::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac){
  bool rv;

  // NOX::Epetra::FiniteDifferenceColoringWithUpdate object
  FiniteDifferenceColoringWithUpdate* testMatrix =
    dynamic_cast<FiniteDifferenceColoringWithUpdate*>(&Jac);
  if (testMatrix == 0) {
    std::cerr << "ERROR: NOX::Epetra::FiniteDifferenceColoringWithUpdate::computeJacobian() - "
     << "Jacobian to evaluate is not a FiniteDifferenceColoringWithUpdate object!"
         << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if(jacobianComputed && use_update){
    rv=differenceProbe(x,*testMatrix->jacobian,*updateColorMap_);
  }
  else{
    rv=differenceProbe(x,*testMatrix->jacobian,*colorMap_);
    jacobianComputed=rv;
  }
  return rv;
}


/****************************************************/
bool FiniteDifferenceColoringWithUpdate::differenceProbe(const Epetra_Vector& x, Epetra_CrsMatrix& jac,const Epetra_MapColoring& colors){

  // Allocate space for perturbation, get column version of x for scaling
  Epetra_Vector xp(x);
  Epetra_Vector *xcol;
  int N=jac.NumMyRows();

  if(jac.ColMap().SameAs(x.Map()))
     xcol=const_cast<Epetra_Vector*>(&x);
  else{
    xcol=new Epetra_Vector(jac.ColMap(),true);//zeros out by default
    xcol->Import(x,*jac.Importer(),InsertAdd);
  }

  // Counters for probing diagnostics
  double tmp,probing_error_lower_bound=0.0,jc_norm=0.0;

  // Grab coloring info (being very careful to ignore color 0)
  int Ncolors=colors.MaxNumColors()+1;
  int num_c0_global,num_c0_local=colors.NumElementsWithColor(0);
  colors.Comm().MaxAll(&num_c0_local,&num_c0_global,1);
  if(num_c0_global>0) Ncolors--;

  if(Ncolors==0) return false;

  // Pointers for Matrix Info
  int entries, *indices;
  double *values;

  // NTS: Fix me
  if ( diffType == Centered ) exit(1);

  double scaleFactor = 1.0;
  if ( diffType == Backward )
    scaleFactor = -1.0;

  // Compute RHS at initial solution
  computeF(x,fo,NOX::Epetra::Interface::Required::FD_Res);

  /* Probing, vector by vector since computeF does not have a MultiVector interface */
  // Assume that anything with Color 0 gets ignored.
  for(int j=1;j<Ncolors;j++){
    xp=x;
    for(int i=0;i<N;i++){
      if(colors[i]==j)
    xp[i] += scaleFactor*(alpha*abs(x[i])+beta);
    }

    computeF(xp, fp, NOX::Epetra::Interface::Required::FD_Res);

    // Do the subtraction to estimate the Jacobian (w/o including step length)
    Jc.Update(1.0, fp, -1.0, fo, 0.0);

    // Relative error in probing
     if(use_probing_diags){
       Jc.Norm2(&tmp);
       jc_norm+=tmp*tmp;
     }

    for(int i=0;i<N;i++){
      // Skip for uncolored row/columns, else update entries
      if(colors[i]==0) continue;

      jac.ExtractMyRowView(i,entries,values,indices);
      for(int k=0;k<jac.NumMyEntries(i);k++){
    if(colors[indices[k]]==j){
      values[k]=Jc[i] / (scaleFactor*(alpha*abs((*xcol)[indices[k]])+beta));
      // If probing diagnostics are on, zero out the entries as they are used
      if(use_probing_diags) Jc[i]=0.0;
      break;// Only one value per row...
    }
      }
    }
    if(use_probing_diags){
      Jc.Norm2(&tmp);
      probing_error_lower_bound+=tmp*tmp;
    }
  }

  // If diagnostics are requested, output Frobenius norm lower bound
  if(use_probing_diags && !x.Comm().MyPID()) printf("Probing Error Lower Bound (Frobenius) abs = %6.4e rel = %6.4e\n",sqrt(probing_error_lower_bound),sqrt(probing_error_lower_bound)/sqrt(jc_norm));

  // Cleanup
  if(!jac.ColMap().SameAs(x.Map()))
    delete xcol;

  return true;
}


