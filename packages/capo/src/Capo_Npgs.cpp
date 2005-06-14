//
// @HEADER
// ***********************************************************************
// 
//                           Capo Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/************************************************************ 
File:      Capo_Npgs.cpp
Purpose:   The Newton--Picard Gauss-Seidel Solver for Capo.
Date:      6-13-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Npgs.hpp"

using namespace CAPO;

//-----------------------------------------------------------------
// Function      : Npgs::Npgs
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Npgs::Npgs(Teuchos::RefCountPtr<Parameter_List> ParamList, \
	   Teuchos::RefCountPtr<Integrator> App_Int, \
	   Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x0, \
	   double lambda0, double T0) 
{
  xcurrent = x0;
  *xfinal = *x0;  //xfinal should have a copy of what x0 is pointing to.

  lambdacurrent = lambda0;
  lambdafinal = lambda0;

  Tcurrent = T0;
  Tfinal = T0;

  iter = 0;

  App_Integrator = App_Int;
  SolveParameters = ParamList;

  // Use the VectorBase xcurrent to create a MultiVector Base.
  Ve = xcurrent->space()->createMembers(30);

};

Npgs::Npgs()
{};

//-----------------------------------------------------------------
// Function      : Npgs::~Npgs
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------

Npgs::~Npgs()
{
};

//-----------------------------------------------------------------
// Function      : Npgs::Initialize
// Purpose       : Create the initial subspace.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Initialize()
{
  for (int i=1;i<((*SolveParameters).get_NumberXtraVecsSubspace())+1;i++)
    {
      Thyra::randomize(0.1,1.1,&*((*Ve).col(i)));
    }

  Orthonormalize(Ve,(*SolveParameters).get_NumberXtraVecsSubspace());

  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_pe;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > We;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re;
  
  Ve_pe = (*Ve).subView(Thyra::Range1D(1,(*SolveParameters).get_NumberXtraVecsSubspace()));
  Se = (*(*Ve_pe).domain()).createMembers((*SolveParameters).get_NumberXtraVecsSubspace());
  Re = (*(*Ve_pe).domain()).createMembers((*SolveParameters).get_NumberXtraVecsSubspace());
  We = (*Ve_pe).clone_mv();

  //Subspace Iterations
  for (int i=1;i<16;i++)
    {
      We = MatVecs(Ve_pe);
      (*Ve_pe).apply(Thyra::TRANS,*We,&*Se);
      SchurDecomp(Se,Re);

      if (i<15)
	{
	  (*We).apply(Thyra::NOTRANS,*Se,&*Ve,1.0);
	  Orthonormalize(Ve,(*SolveParameters).get_NumberXtraVecsSubspace());
	}
    }
  Orthonormalize(Ve,(*SolveParameters).get_NumberXtraVecsSubspace());

  cout << "Npgs::Initialize We have created the initial subspace. " << endl;

}
//-----------------------------------------------------------------
// Function      : Npgs::Orthonormalize
// Purpose       : Orthonormalize the Vectors of a MultiVectorBase.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Orthonormalize(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Basis, int Number_of_Columns)
{
  double temp = 0;
  for(int i=1;i<Number_of_Columns+1;i++)
    {
      for(int j=1;j<i;j++)
	{
	  temp = Thyra::dot(*((*Basis).col(i)),*((*Basis).col(j)));
	  Thyra::Vp_StV( &*((*Basis).col(i)),-temp, *((*Basis).col(j)) );
	}
      temp = Thyra::norm( *((*Basis).col(i)) );
      if (temp > 1.0e-12)
	Thyra::Vt_S( &*((*Basis).col(i)), 1/temp );
      else
	cout <<"WARNING Npgs::Orthonormalize: norm=0; No scaling done" << endl;
    }
}
//-----------------------------------------------------------------
// Function      : Npgs::Finish
// Purpose       : Print the Basis Matrix
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Finish()
{
  if ((*SolveParameters).get_printproc() >0)
    cout <<"\n Writing columns of Ve:" << endl;
  Print(Ve);
}
//-----------------------------------------------------------------
// Function      : Npgs::InnerFunctions
// Purpose       : Npgs does not require anything to be done 
//                 between inner iterations.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::InnerFunctions()
{
}
//-----------------------------------------------------------------
// Function      : Npgs::Predictor
// Purpose       : Currently increments lambda, and sets the new
//                 solution guesses to the current solution.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Predictor()
{
  lambdacurrent = lambdafinal+(*SolveParameters).get_lambda_stepsize();
  Tcurrent = Tfinal;
  *xcurrent = *xfinal;
}
//-----------------------------------------------------------------
// Function      : Npgs::Get_Tfinal
// Purpose       : Returns Tfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& Npgs::Get_Tfinal()
{
  return Tfinal;
}
//-----------------------------------------------------------------
// Function      : Npgs::Get_lambdafinal
// Purpose       : Returns lambdafinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& Npgs::Get_lambdafinal()
{
  return lambdafinal;
}
//-----------------------------------------------------------------
// Function      : Npgs::Get_xfinal
// Purpose       : Returns xfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Npgs::Get_xfinal()
{
  return xfinal;
}
//-----------------------------------------------------------------
// Function      : Npgs::MatVec
// Purpose       : returns M*y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > Npgs::MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > y)
{
  double delta = 1.0e-5; /*Finite Difference constant */
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > upert;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phiupert;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phiu;

  Thyra::assign(&*upert, 0.0);
  Thyra::Vp_StV(&*upert,delta,*y);
  Thyra::Vp_StV(&*upert,1.0,*xcurrent);

  (*App_Integrator).Integrate(xcurrent,phiu,Tcurrent,lambdacurrent);
  (*App_Integrator).Integrate(upert,phiupert,Tcurrent,lambdacurrent);

  Thyra::Vp_StV(&*phiupert,-1.0,*phiu);
  Thyra::Vt_S(&*phiupert,1/delta);

  return phiupert;
}
//-----------------------------------------------------------------
// Function      : Npgs::MatVecs
// Purpose       : returns M*Y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Npgs::MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Y)
{
  int dimension = (*(*Y).domain()).dim();
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > result;
  result = (*(*xcurrent).space()).createMembers(dimension);
 
  for (int i=1;i<dimension+1;i++)
    {
      (*result).col(i)=MatVec((*Y).col(i));
    }
  return result;
}
//-----------------------------------------------------------------
// Function      : Npgs::IterationFailure
// Purpose       : If the inner iteration function fails to find a
//                 solution, then the code must exit.  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::IterationFailure()
{
  cout << "Npgs::IterationFailure Failed to find a solution." << endl;
  cout << "I should find a better way of doing this..." << endl;
  abort();
}
//-----------------------------------------------------------------
// Function      : Npgs::InnerIteration
// Purpose       : The Npgs algorithm  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
bool Npgs::InnerIteration()
{

  cout << "Inside Inner Iteration loop!" << endl;
  return true;

}
//-----------------------------------------------------------------
// Function      : Npgs::SchurDecomposition
// Purpose       : One Multivector comes in, Schur Decomposition
//                 is performed on it, and the vectors are returned
//                 along with a mutlivector containing the scaling factors.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::SchurDecomp(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se,Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re)
{
  //Fill all this in later!!
  cout << "The Schur decomposition function has not yet been implemented!" << endl;
}
//-----------------------------------------------------------------
// Function      : Npgs::Print
// Purpose       : Print a MultiVector
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Print(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Printme)
{
  cout << "I have not yet implemented the Print function!" << endl;
}
