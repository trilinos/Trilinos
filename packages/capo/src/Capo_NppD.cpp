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
File:      Capo_NppD.cpp
Purpose:   The Newton--Picard Gauss-Seidel Solver for Capo.
Date:      7-29-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_NppD.hpp"

using namespace CAPO;

int Select(double *x, double *y);
int NegSelect(double *x, double *y);

//-----------------------------------------------------------------
// Function      : NppD::NppD
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
NppD::NppD(const Teuchos::RefCountPtr<Parameter_List>& ParamList, 
	   const Teuchos::RefCountPtr<Integrator>& App_Int, 
	   const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x0, 
	   double lambda0, double T0) 
{
  //Note, if we are using this algorithm, both continuation
  // and periodic should be turned on.  I've removed
  // any checks from within the code and assume
  // the variables are set to true. 

  xprevious = createMember(x0->space());
  xcurrent = createMember(x0->space());
  xfinal = createMember(x0->space());
  xinit = createMember(x0->space());
  finit = createMember(x0->space());
  xstep = createMember(x0->space());

  Thyra::assign(&*xprevious, *x0); 
  Thyra::assign(&*xcurrent, *x0); 
  Thyra::assign(&*xfinal, *x0);
  Thyra::assign(&*xstep, 0.0);

  lambdaprevious = lambda0;
  lambdacurrent = lambda0;
  lambdafinal = lambda0;
  lambdastep = 0.0;

  Tprevious = T0;
  Tcurrent = T0;
  Tfinal = T0;
  Tstep = 0.0;

  /** For period doubling bifurcation **/
  dvector = createMember(x0->space());
  Createvectord();
  nuprevious = createMember(x0->space());
  nucurrent = createMember(x0->space());
  nufinal = createMember(x0->space());
  nustep = createMember(x0->space());
  Thyra::assign(&*nustep, 0.0);
  Thyra::assign(&*nuprevious, *dvector); 
  Thyra::assign(&*nucurrent, *dvector); 
  Thyra::assign(&*nufinal, *dvector);

  App_Integrator = App_Int;
  SolveParameters = ParamList;

  Thyra::assign(&*xinit,*xcurrent);
  App_Integrator->Integrate(xcurrent,finit,10e-8,lambdacurrent);
  Thyra::Vp_StV(&*finit,-1.0,*xcurrent);
  Thyra::Vt_S(&*finit,1.0/10e-8);

  iter = 0;
  Unstable_Basis_Size = 0;


  First_Continuation_Step = true;

  // Use the VectorBase xcurrent to create a MultiVector Base.
  Ve = Thyra::createMembers(xcurrent->space(),30);

}

//-----------------------------------------------------------------
// Function      : NppD::Createvectord
// Purpose       : Create the initial subspace.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
void NppD::Createvectord()
{
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_pe;
  Ve_pe = createMembers(xcurrent->space(),SolveParameters->get_NumberXtraVecsSubspace());

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVector;
  TempVector = createMember(xcurrent->space());

  for (int i=1;i<SolveParameters->get_NumberXtraVecsSubspace()+1;i++)
    {
      
      Thyra::assign(&*TempVector,*Ve_pe->col(i));
      Thyra::randomize(0.1,1.1,&*TempVector);
      Thyra::assign(&*Ve_pe->col(i),*TempVector);
      
    }
  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());

  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > We;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re;
  
  Se = createMembers(Ve_pe->domain(),SolveParameters->get_NumberXtraVecsSubspace());
  Re = createMembers(Ve_pe->domain(),SolveParameters->get_NumberXtraVecsSubspace());
  We = createMembers(xcurrent->space(),SolveParameters->get_NumberXtraVecsSubspace());

  //Subspace Iterations
  for (int i=0;i<15;i++)
    {
      We = MatVecs(Ve_pe);
      Thyra::apply(*Ve_pe,Thyra::TRANS,*We,&*Se);
      SchurDecomp(Se,Re,1);
      if (i<14)
	{
	  Thyra::apply(*We,Thyra::NOTRANS,*Se,&*Ve_pe);
	  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());
	}
    }
  
  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());

  Thyra::assign( &*dvector,*(Ve_pe->col(1)) );
  // The Schur decomposition may have altered the Unstable basis size
  // so reset it to zero.
  Unstable_Basis_Size = 0;
  cout << "Vector d/v created." << endl;

}

//-----------------------------------------------------------------
// Function      : NppD::Initialize
// Purpose       : Create the initial subspace.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void NppD::Initialize()
{
  unsigned int seedtime = time(NULL);
  Teuchos::ScalarTraits<Scalar>::seedrandom(seedtime);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVector;
  TempVector = createMember(xcurrent->space());

  for (int i=1;i<30+1;i++)
    {
      
      Thyra::assign(&*TempVector,*Ve->col(i));
      Thyra::randomize(0.1,1.1,&*TempVector);
      Thyra::assign(&*Ve->col(i),*TempVector);
      
    }

  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_pe;

  Ve_pe = createMembers(xcurrent->space(),SolveParameters->get_NumberXtraVecsSubspace());
  Ve_pe = Ve->subView(Thyra::Range1D(1,SolveParameters->get_NumberXtraVecsSubspace()));

  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());

  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > We;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re;
  
  Se = createMembers(Ve_pe->domain(),SolveParameters->get_NumberXtraVecsSubspace());
  Re = createMembers(Ve_pe->domain(),SolveParameters->get_NumberXtraVecsSubspace());
  We = createMembers(xcurrent->space(),SolveParameters->get_NumberXtraVecsSubspace());

  //Subspace Iterations
  for (int i=0;i<15;i++)
    {
      We = MatVecs(Ve_pe);
      Thyra::apply(*Ve_pe,Thyra::TRANS,*We,&*Se);
      SchurDecomp(Se,Re,0);
      if (i<14)
	{
	  Thyra::apply(*We,Thyra::NOTRANS,*Se,&*Ve_pe);
	  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());
	}
    }
  
  Orthonormalize(Ve_pe,SolveParameters->get_NumberXtraVecsSubspace());
  // The Schur decomposition may have altered the Unstable basis size
  // so reset it to zero.
  Unstable_Basis_Size = 0;
  cout << "Initial subspace created. " << endl;

}


//-----------------------------------------------------------------
// Function      : NppD::Orthonormalize
// Purpose       : Orthonormalize the Vectors of a MultiVectorBase.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/23/05
//------------------------------------------------------------------
void NppD::Orthonormalize(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Basis, int Number_of_Columns)
{
  double temp = 0;
  temp = sqrt( Thyra::dot( *(Basis->col(1)),*(Basis->col(1))));
  if (temp > 1.0e-12)
    Thyra::Vt_S( &*(Basis->col(1)), 1.0/temp );
  else
    cout <<"WARNING NppD::Orthonormalize: norm=0; No scaling done" << endl;


  for(int i=2;i<Number_of_Columns+1;i++)
    {
      for(int j=1;j<i;j++)
	{
	  temp = ( Thyra::dot(*(Basis->col(i)),*(Basis->col(j))) );
	  Thyra::Vp_StV( &*(Basis->col(i)), -temp, *(Basis->col(j)) );
	}
      temp = ( Thyra::dot( *(Basis->col(i)),*(Basis->col(i)) ) );
      if (temp > 1.0e-12)
	Thyra::Vt_S( &*(Basis->col(i)), 1.0/sqrt(temp) );
      else
	cout <<"WARNING NppD::Orthonormalize: norm=0; No scaling done" << endl;
    }
}


//-----------------------------------------------------------------
// Function      : NppD::Finish
// Purpose       : Print the Basis Matrix
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void NppD::Finish()
{
  if ((*SolveParameters).get_printproc() >0)
    cout << "\n NPGS finished " << endl;
    //cout <<"\n Writing columns of Ve:" << endl;
  //Print(Ve);
}
//-----------------------------------------------------------------
// Function      : NppD::InnerFunctions
// Purpose       : NppD does not require anything to be done 
//                 between inner iterations.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void NppD::InnerFunctions()
{
}
//-----------------------------------------------------------------
// Function      : NppD::Predictor
// Purpose       : Currently increments lambda, and sets the new
//                 solution guesses to the current solution.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/20/05
//------------------------------------------------------------------
void NppD::Predictor()
{

  // Eventually put into parameter list...
  double lambdamax = SolveParameters->get_lambda_max();
  double lambdamin = SolveParameters->get_lambda_min();

  Tstep = Tfinal-Tprevious;
  lambdastep = lambdafinal-lambdaprevious;
  Thyra::assign(&*xstep,*xfinal);
  Thyra::Vp_StV(&*xstep,-1.0,*xprevious);
  Thyra::assign(&*nustep,*nufinal);
  Thyra::Vp_StV(&*nustep,-1.0,*nuprevious);

  Tprevious = Tfinal;
  lambdaprevious = lambdafinal;
  Thyra::assign(&*xprevious,*xfinal);
  Thyra::assign(&*nuprevious,*nufinal);
  if ( SolveParameters->Arc_Length() )
    {
      if ( iter < SolveParameters->get_lambda_extend_tol() )
	lambdastep = 2*lambdastep;
      if (lambdastep > lambdamax)
	lambdastep = lambdamax;
      if (lambdastep < lambdamin)
	lambdastep = lambdamin;
    }
  if ( First_Continuation_Step )
    {
      First_Continuation_Step = false;
      Tcurrent = Tprevious;
      lambdacurrent = lambdaprevious+(*SolveParameters).get_lambda_stepsize();
      Thyra::assign(&*xcurrent,*xprevious);
      Thyra::assign(&*nucurrent,*nuprevious);
    }
  else
    {
      Tcurrent = Tprevious + Tstep;
      lambdacurrent = lambdaprevious + lambdastep;

      Thyra::assign(&*xcurrent,*xprevious);
      Thyra::Vp_StV(&*xcurrent,1.0,*xstep);
      Thyra::assign(&*nucurrent,*nuprevious);
      Thyra::Vp_StV(&*nucurrent,1.0,*nustep);
    }

}
//-----------------------------------------------------------------
// Function      : NppD::Get_Tfinal
// Purpose       : Returns Tfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& NppD::Get_Tfinal()
{
  return Tfinal;
}


//-----------------------------------------------------------------
// Function      : NppD::Get_lambdafinal
// Purpose       : Returns lambdafinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& NppD::Get_lambdafinal()
{
  return lambdafinal;
}


//-----------------------------------------------------------------
// Function      : NppD::Get_xfinal
// Purpose       : Returns xfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& NppD::Get_xfinal()
{
  return xfinal;
}


//-----------------------------------------------------------------
// Function      : NppD::MatVec
// Purpose       : returns M*y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > NppD::MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y)
{
  double delta = 1.0e-5; /*Finite Difference constant */
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > upert;
  upert = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phiupert;
  phiupert = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phiu;
  phiu = createMember(xcurrent->space());

  Thyra::assign(&*upert, 0.0);
  Thyra::Vp_StV(&*upert,delta,*y);
  Thyra::Vp_StV(&*upert,1.0,*xfinal);

  App_Integrator->Integrate(xfinal,phiu,Tfinal,lambdafinal);
  App_Integrator->Integrate(upert,phiupert,Tfinal,lambdafinal);

  Thyra::Vp_StV(&*phiupert,-1.0,*phiu);
  Thyra::Vt_S(&*phiupert,1.0/delta);

  return phiupert;
}


//-----------------------------------------------------------------
// Function      : NppD::MatVecs
// Purpose       : returns M*Y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > NppD::MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Y)
{
  int dimension = Y->domain()->dim();
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > result;
  result = createMembers(xcurrent->space(),dimension);
  for (int i=1;i<dimension+1;i++)
    {
      Thyra::assign(&*(result->col(i)),*(MatVec(Y->col(i))));
    }
  return result;
}


//-----------------------------------------------------------------
// Function      : NppD::IterationFailure
// Purpose       : If the inner iteration function fails to find a
//                 solution, then the code must exit.  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void NppD::IterationFailure()
{
  cout << "NppD::IterationFailure Failed to find a solution." << endl;
  lambdastep = .5*lambdastep;
  Tstep = .5*Tstep;
  Thyra::Vt_S(&*xstep,.5);
  Thyra::Vt_S(&*nustep,.5);

  lambdacurrent = lambdaprevious + lambdastep;
  Tcurrent = Tprevious + Tstep;
  Thyra::assign(&*xcurrent,*xprevious);
  Thyra::Vp_StV(&*xcurrent,1.0,*xstep);
  Thyra::assign(&*nucurrent,*nuprevious);
  Thyra::Vp_StV(&*nucurrent,1.0,*nustep);

  cout << "Cutting the lambdastep in half." << endl;
  //abort();
}


//-----------------------------------------------------------------
// Function      : NppD::InnerIteration
// Purpose       : The NppD algorithm  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool NppD::InnerIteration()
{

  /** Declarations **/
  bool converged=false;

  double eps = 10e-8;
  double *deltaT1;
  deltaT1 = new double;
  double *deltaT2;
  deltaT2 = new double;
  double *deltalambda;
  deltalambda = new double;

  int Subspace_Size = 0;


  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > v;
  v = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > r;
  r = createMember(xcurrent->space());


  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xq1;
  xq1 = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xq2;
  xq2 = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > nuq1;
  nuq1 = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > nuq2;
  nuq2 = createMember(xcurrent->space());

  App_Integrator->Integrate(xcurrent,v,Tcurrent,lambdacurrent);
  Thyra::assign(&*r,*v); 
  Thyra::Vp_StV(&*r,-1.0,*xcurrent);

  // xfinal, Tfinal, and lambdafinal change... but start off at
  // same values as currents...
  Thyra::assign(&*xfinal,*xcurrent); 
  Tfinal = Tcurrent;
  lambdafinal = lambdacurrent;

  converged = Converged(xcurrent,v);
  iter = 0;
  while ((iter<SolveParameters->get_MaxInnerIts()) && (!converged))
    {
      //Declarations
      Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_pe;
      Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se;
      Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > We;
      Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re;
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVector;
      TempVector = createMember(xcurrent->space());

      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dphidlambda;
      dphidlambda = createMember(xcurrent->space());
      dphi_dlambda(dphidlambda);

      // Setup
      Subspace_Size = SolveParameters->get_NumberXtraVecsSubspace()+Unstable_Basis_Size;

      Ve_pe = createMembers(xcurrent->space(),Subspace_Size);
      Ve_pe = Ve->subView(Thyra::Range1D(1,Subspace_Size));
      Se = createMembers(Ve_pe->domain(),Subspace_Size);
      Re = createMembers(Ve_pe->domain(),Subspace_Size);
      We = createMembers(xcurrent->space(),Subspace_Size);
      

      // Algorithm
      SubspaceIterations(Se,We,Re);
      if (Unstable_Basis_Size > 0) // Need to include a Newton-Step.
	{
	  // Declarations (these change size with every go-around)
	  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Vp;
	  Vp = createMembers(xcurrent->space(),Unstable_Basis_Size);

	  ComputeVp(Se,Vp);

	  /** Calculate xq's **/
	  Calculatedq(Vp,xq1,r);
	  Calculatedq(Vp,xq2,dphidlambda);

	  /** Calculate xp's **/
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xp1;
	  xp1 = createMember(Vp->domain());
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xp2;
	  xp2 = createMember(Vp->domain());
	  Calculatedp(Vp,xq1,xq2,xp1,xp2,dphidlambda,Re,v,finit,r,*deltaT1,*deltaT2);

	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xi1;
	  xi1 = createMember(xcurrent->space());
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xi2;
	  xi2 = createMember(xcurrent->space());
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > Mpartialvec;
	  Mpartialvec = createMember(xcurrent->space());

	  MT(Mpartialvec,nufinal,v);

	  /** Build the xi's **/
	  xi1 = MatVec(nufinal);
	  Thyra::Vp_StV(&*xi1,*deltaT1,*Mpartialvec);
	  Thyra::Vp_StV(&*xi1,1.0,*nufinal);

	  ML(xi2,nufinal);
	  Thyra::Vp_StV(&*xi2,*deltaT2,*Mpartialvec);

	  Thyra::apply(*Vp,Thyra::NOTRANS,*xp1,&*TempVector);
	  Thyra::Vp_StV(&*TempVector,1.0,*xq1);
	  MX(Mpartialvec,nufinal,TempVector);
	  Thyra::Vp_StV(&*xi1,1.0,*Mpartialvec);

	  Thyra::apply(*Vp,Thyra::NOTRANS,*xp2,&*TempVector);
	  Thyra::Vp_StV(&*TempVector,1.0,*xq2);
	  MX(Mpartialvec,nufinal,TempVector);
	  Thyra::Vp_StV(&*xi2,1.0,*Mpartialvec);

	  /** Calculate nuq's **/
	  Calculatedq(Vp,nuq1,xi1);
	  Calculatedq(Vp,nuq2,xi2);
	  
	  /** Calculate nup's **/
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > nup1;
	  nup1 = createMember(Vp->domain());
	  Calculatenup(Vp,nuq1,nuq2,nup1,dphidlambda,Re,xi1,xi2,*deltalambda);

	  /** Update dt,dlambda **/
	  lambdafinal+= *deltalambda;
	  Tfinal+=(*deltaT1+ *deltalambda* (*deltaT2) );
	  /** Update xfinal **/
	  Thyra::Vp_StV(&*xfinal,1.0,*xq1);
	  Thyra::apply(*Vp,Thyra::NOTRANS,*xp1,&*TempVector);
	  Thyra::Vp_StV(&*xfinal,1.0,*TempVector);
	  Thyra::apply(*Vp,Thyra::NOTRANS,*xp2,&*TempVector);
	  Thyra::Vp_StV(&*TempVector,1.0,*xq2);
	  Thyra::Vp_StV(&*xfinal,*deltalambda,*TempVector);
	  /** Update nufinal **/
	  Thyra::Vp_StV(&*nufinal,1.0,*nuq1);
	  Thyra::Vp_StV(&*nufinal,*deltalambda,*nuq2);
	  Thyra::apply(*Vp,Thyra::NOTRANS,*nup1,&*TempVector);
	  Thyra::Vp_StV(&*nufinal,1.0,*TempVector);

	  UpdateVe(We,Se);
	}
      else
	{
	  if ( SolveParameters->Periodic() )
	    {
	      cout << "Error: If the Problem is periodic, there should" << endl;
	      cout << "       always be at least one floquet multiplier " << endl;
	      cout << "       of value 1>rho!" << endl;
	      cout << "The code will assume deltaT=0, deltalambda=0 and continue. " << endl;
	      
	      TempVector = MatVec(r);
	      Thyra::Vp_StV(&*TempVector,1.0,*r);
	      TempVector = MatVec(TempVector);
	      Thyra::Vp_StV(&*TempVector,1.0,*r);
	      Thyra::Vp_StV(&*xfinal,1.0,*TempVector);

	    }
	  
	}

      Orthonormalize(Ve,Unstable_Basis_Size+SolveParameters->get_NumberXtraVecsSubspace());
      App_Integrator->Integrate(xfinal,v,Tfinal,lambdafinal);
      Thyra::assign(&*r,*v); 
      Thyra::Vp_StV(&*r,-1.0,*xfinal);

      converged = Converged(xfinal,v);


      iter++;
    }//end while

  delete deltaT1;
  delete deltaT2;
  delete deltalambda;

  return converged;

}


//-----------------------------------------------------------------
// Function      : NppD::SchurDecomposition
// Purpose       : One Multivector comes in, Schur Decomposition
//                 is performed on it, and the vectors are returned
//                 along with a multivector containing the scaling factors.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
void NppD::SchurDecomp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re, int flag)
{
  /* This function performs a Schur Decomposition.
     On return Mat contains the Upper Triangular matrix with 
     eigenvalues on the diagonal, and V contains the basis
     vectors.  This function also requires the existence of
     function SELECT which is called in when ordering the 
     eigenvalues along the diagonal.
  */
  int m = Se->domain()->dim();
  
  double *se = new double[m*m];
  double *re = new double[m*m];

  RTOpPack::SubMultiVectorT<Scalar> sub_mv;
  Se->getSubMultiVector( Thyra::Range1D(1,m) , Thyra::Range1D(1,m),&sub_mv);
  for (int j=0;j<sub_mv.numSubCols();j++)
    {
      for (int i=0;i<sub_mv.subDim();i++)
	{
	  se[i+j*m] = sub_mv.values()[i+j*sub_mv.leadingDim()];
	  re[i+j*m] = 0;
	}
    }
  Se->freeSubMultiVector(&sub_mv);
  int lwork = 10*m+10;
  char *cn="V";
  char *cs="S"; //Sort the eigenvalues using the SELECT function.
  int sdim=0;
  double *wr = new double[m]; //These contain the real and imag. parts
  double *wi = new double[m]; // of the eigenvalues.
  double *work = new double [lwork];
  int *bwork = new int[m];
  int info=0;
  int LDA = m;
  int (*pt2Select)(double*, double*);


  if (flag == 0)
    {
      pt2Select = &Select;
    }
  if (flag == 1)
    {
      pt2Select = &NegSelect;
    }

  dgees_(cn, cs, pt2Select, &m, se, &LDA, &sdim, 
	 wr, wi, re, &LDA, work, &lwork, bwork, &info);

  /* On output the array re contains the schur vectors which I want
     in the matrix Se.
  */

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > tempcol1;
  tempcol1 = createMember(Se->domain());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > tempcol2;
  tempcol2 = createMember(Re->domain());

  for (int i=0;i<m;i++)
    {
      for (int j=0;j<m;j++)
	{
	  Thyra::set_ele(j+1,se[j+i*m],&*tempcol1);
	  Thyra::set_ele(j+1,re[j+i*m],&*tempcol2);
	}
      Thyra::assign(&*(Se->col(i+1)),*tempcol2);
      Thyra::assign(&*(Re->col(i+1)),*tempcol1);
    }

  /* sdim tells us the number of eigenvalues that met the
     criteria of the selection function, in the case of 
     the Newton-Picard algorithm, its the number of Floquet
     multipliers outside of disk of specified radius. 
  */
  Unstable_Basis_Size = sdim;

  delete [] work;
  delete [] bwork;
  delete [] wr;
  delete [] wi;
  delete [] se;
  delete [] re;
  if (info != 0)
    cout << "There was an error returned from the Schur decompostion." << endl;

}


//-----------------------------------------------------------------
// Function      : NppD::Print
// Purpose       : Print a MultiVector
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void NppD::Print(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Printme)
{
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TmpColumn;
  TmpColumn = createMember(xcurrent->space());

  int Column_Number = Printme->domain()->dim();
  int Row_Number = Printme->range()->dim();
  for (int i=1;i<Column_Number+1;i++)
    {
      Thyra::assign(&*TmpColumn,*(Printme->col(i)));
      cout << "Column " << i << endl;
      for (int j=1;j<Row_Number+1;j++)
	{
	  cout << Thyra::get_ele(*TmpColumn,j) << endl;
	}
    }

}


//-----------------------------------------------------------------
// Function      : NppD::Select
// Purpose       : A necessary function for the Schur Decomposition
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
int Select(double *x, double *y)
{
  /* This function is only called by the ScurDecomp function.  It 
     is used as the selection criterion when sorting eigenvalues.  
     It is supposed to return a true or false based on whether or
     not the eigenvalue x+iy meets a specified selection criterion.
     Because it will be called by a fortran code, it must return 
     an int rather than a bool.  Thus, false=0, true=1.  For now
     I'll just use a straight magnitude criterion.
  */

  double z;
  z = sqrt((*x)*(*x)+(*y)*(*y));
  if (z>.5)
    return 1;
  else
    return 0;
}


//-----------------------------------------------------------------
// Function      : NppD::NegSelect
// Purpose       : A necessary function for the Schur Decomposition
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
int NegSelect(double *x, double *y)
{
  /* This function is only called by the ScurDecomp function.  It 
     is used as the selection criterion when sorting eigenvalues.  
     It is supposed to return a true or false based on whether or
     not the eigenvalue x+iy meets a specified selection criterion.
     Because it will be called by a fortran code, it must return 
     an int rather than a bool.  Thus, false=0, true=1.  For now
     I'll just use a straight magnitude criterion.
  */

  double z;
  z = sqrt((*x)*(*x)+(*y)*(*y));
  if ((*x > -1.01) && (*x < -.99))
    return 1;
  else
    return 0;
}


//-----------------------------------------------------------------
// Function      : NppD::Converged
// Purpose       : Checks to see if the norm(x-y)<tol
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool NppD::Converged(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,
	       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y)
{
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > differencevector;
  differencevector = createMember(x->space());
  int n = y->space()->dim();

  Thyra::assign(&*differencevector,*x);
  Thyra::Vp_StV(&*differencevector,-1.0,*y);
  double diffnorm = (sqrt( Thyra::dot(*differencevector,*differencevector) ) )/ (double) n;
  cout << "The norm of the residual is " << diffnorm << endl;
  if (diffnorm < SolveParameters->get_tol())
    return true;
  else
    return false;
}

//-----------------------------------------------------------------
// Function      : NppD::SubspaceIterations
// Purpose       : See Lust et.al.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void NppD::SubspaceIterations(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We, const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re)
{
  int Subspace_Size = SolveParameters->get_NumberXtraVecsSubspace()+Unstable_Basis_Size;
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_pe;

  Ve_pe = createMembers(xcurrent->space(),Subspace_Size);
  Ve_pe = Ve->subView(Thyra::Range1D(1,Subspace_Size));
  Orthonormalize(Ve_pe,Subspace_Size);

  for (int i=0;i<SolveParameters->get_SubspaceIterations();i++)
    {
      We = (MatVecs(Ve_pe));
      Thyra::apply(*Ve_pe,Thyra::TRANS,*We,&*Se);
      SchurDecomp(Se,Re,0);
      if (i<(SolveParameters->get_SubspaceIterations()-1))
	{
	  Thyra::apply(*We,Thyra::NOTRANS,*Se,&*Ve_pe);
	  Orthonormalize(Ve_pe,Subspace_Size);
	}
    }
  Orthonormalize(Ve_pe,Subspace_Size);
}

//-----------------------------------------------------------------
// Function      : NppD::Calculatedq
// Purpose       : Perform the Picard iterations
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void NppD::Calculatedq(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp
		       ,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r)
{
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > q;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVec1;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVec2;

  q = createMember(r->space());
  TempVec1 = createMember(Vp->domain());
  TempVec2 = createMember(r->space());

  Thyra::apply(*Vp,Thyra::TRANS,*r,&*TempVec1);
  Thyra::apply(*Vp,Thyra::NOTRANS,*TempVec1,&*q);

  Thyra::Vt_S(&*q,-1.0);
  Thyra::Vp_StV(&*q,1.0,*r);

  TempVec2 = MatVec(q);
  Thyra::Vp_StV(&*TempVec2,1.0,*r);
  Thyra::apply(*Vp,Thyra::TRANS,*TempVec2,&*TempVec1);
  Thyra::apply(*Vp,Thyra::NOTRANS,*TempVec1,&*q);

  Thyra::assign(&*dq,*TempVec2);
  Thyra::Vp_StV(&*dq,-1.0,*q);

  Thyra::assign(&*q,*dq);
  TempVec2 = MatVec(q);
  Thyra::Vp_StV(&*TempVec2,1.0,*r);
  Thyra::apply(*Vp,Thyra::TRANS,*TempVec2,&*TempVec1);
  Thyra::apply(*Vp,Thyra::NOTRANS,*TempVec1,&*q);

  Thyra::assign(&*dq,*TempVec2);
  Thyra::Vp_StV(&*dq,-1.0,*q);

}


//-----------------------------------------------------------------
// Function      : NppD::ComputeVp
// Purpose       : Extract Vp from Ve and multiply it by Se.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool NppD::ComputeVp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,
		     const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp)
{
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve_p;
  Ve_p = createMembers(xcurrent->space(),Unstable_Basis_Size);
  Ve_p = Ve->subView(Thyra::Range1D(1,Unstable_Basis_Size));
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Yp;
  Yp = createMembers(Vp->domain(),Unstable_Basis_Size);
 
  RTOpPack::SubMultiVectorT<Scalar> sub_se;
  RTOpPack::MutableSubMultiVectorT<Scalar> sub_yp;
  
  Se->getSubMultiVector(Thyra::Range1D(1,Unstable_Basis_Size),Thyra::Range1D(1,Unstable_Basis_Size),&sub_se);
  Yp->getSubMultiVector(Thyra::Range1D(1,Unstable_Basis_Size),Thyra::Range1D(1,Unstable_Basis_Size),&sub_yp);
  
  for (int j=0;j<sub_se.numSubCols();j++)
    {
      for (int i=0;i<sub_se.subDim();i++)
	{
	  sub_yp.values()[i+j*sub_yp.leadingDim()] = 
	    sub_se.values()[i+j*sub_se.leadingDim()];
	}
    }
  Se->freeSubMultiVector(&sub_se);
  Yp->commitSubMultiVector(&sub_yp);

  Thyra::apply(*Ve_p,Thyra::NOTRANS,*Yp,&*Vp); //Vp=Ve*Se[1..p,1..p]
  return true;
}


//-----------------------------------------------------------------
// Function      : NppD::UpdateVe
// Purpose       : Update the Basis
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool NppD::UpdateVe(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We,
		    const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se)
{
  int Size = We->domain()->dim();
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > TempMV;
  TempMV = createMembers(xcurrent->space(),Size);
  TempMV = Ve->subView(Thyra::Range1D(1,Size));
  Thyra::apply(*We,Thyra::NOTRANS,*Se,&*TempMV);
  return true;

}
//-----------------------------------------------------------------
// Function      : NppD::Calcuatedp
// Purpose       : Find the xp1,xp2 vectors and deltaT1,deltaT2
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/28/05
//------------------------------------------------------------------
bool NppD::Calculatedp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xq1,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xq2,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xp1,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xp2,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dphidlambda,
		       const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
		       double& deltaT1,double& deltaT2)
{
  double eps = 10e-5;

  // Dimension of the linear system depends on whether we are looking 
  // for a periodic solution.
  int LS_size = Unstable_Basis_Size+1;
    
  // First need a square MultiVector...
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > TempMV;
  TempMV = createMembers(xcurrent->space(),LS_size);
  double *LHS;
  double *RHS;
  LHS = new double[(LS_size)*(LS_size)];
  RHS = new double[2*LS_size];

  /*Upper left corner of the lhs */
  for (int j=0;j<Unstable_Basis_Size;j++)
    for (int i=0;i<Unstable_Basis_Size;i++)
      {
	LHS[i+j*(LS_size)]=Thyra::get_ele(*Re->col(j+1),i+1);
	if (i==j)
	  LHS[i+j*(LS_size)]+=-1.0;
      }

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp;
  temp = createMember(xcurrent->space());

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dphidt;
  dphidt = createMember(xcurrent->space());
  dphi_dt(dphidt);

  //Rightmost column of LHS.
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > rightcol;
  rightcol = createMember(TempMV->domain());
  Thyra::apply(*Vp,Thyra::TRANS,*dphidt,&*rightcol);

  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      LHS[i+(Unstable_Basis_Size)*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
    }
  
  //add a row.
  Thyra::apply(*Vp,Thyra::TRANS,*finit,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      LHS[Unstable_Basis_Size+i*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
    }
  
  // and the corner
  LHS[(LS_size)*(LS_size)-1]=0.0;
  
  // Now the Right Hand Side Vector
  temp = MatVec(xq1);
  Thyra::Vp_StV(&*temp,1.0,*r);
  Thyra::apply(*Vp,Thyra::TRANS,*temp,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i]=-Thyra::get_ele(*rightcol,i+1);
    }
  temp = MatVec(xq2);
  Thyra::Vp_StV(&*temp,1.0,*dphidlambda);
  Thyra::apply(*Vp,Thyra::TRANS,*temp,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i+LS_size]=-Thyra::get_ele(*rightcol,i+1);
    }

  Thyra::assign(&*temp,*xfinal);
  Thyra::Vp_StV(&*temp,-1.0,*xinit);
  Thyra::Vp_StV(&*temp,1.0,*xq1);
  RHS[Unstable_Basis_Size]=-Thyra::dot(*finit,*temp);

  Thyra::assign(&*temp,*xq2);
  RHS[2*LS_size-1]=-Thyra::dot(*finit,*temp);

  Solve_Linear(LHS,RHS, false ,LS_size,2);
  for (int j=0;j<Unstable_Basis_Size;j++)
    {
      Thyra::set_ele(j+1,RHS[j],&*xp1);
      Thyra::set_ele(j+1,RHS[j+LS_size],&*xp2);
    }

  deltaT1 = RHS[Unstable_Basis_Size];
  deltaT2 = RHS[2*LS_size+1];

  delete [] LHS;
  delete [] RHS;

  return true;
}


//-----------------------------------------------------------------
// Function      : NppD::dphi_dt
// Purpose       : Use finite differences to calculate f(phi).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool NppD::dphi_dt(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f)
{
  double delta = 1.0e-5; /*Finite Difference constant*/
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phi;
  phi = createMember(f->space());

  App_Integrator->Integrate(xfinal,f,Tfinal+delta,lambdafinal);
  App_Integrator->Integrate(xfinal,phi,Tfinal,lambdafinal);
  Thyra::Vp_StV(&*f,-1.0,*phi);
  Thyra::Vt_S(&*f,1.0/delta);

  return true;
}

//-----------------------------------------------------------------
// Function      : NppD::feval
// Purpose       : Use finite differences to calculate f(x).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
bool NppD::feval(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f,
		 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& u)
{
  double delta = 1.0e-5; /*Finite Difference constant*/
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp;
  temp = createMember(f->space());

  App_Integrator->Integrate(u,f,delta,lambdafinal);
  Thyra::Vp_StV(&*f,-1.0,*u);
  Thyra::Vt_S(&*f,1.0/delta);

  return true;
}
//-----------------------------------------------------------------
// Function      : NppD::Solve_Linear
// Purpose       : Call LAPACK routine to solve a linear system
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/20/05
//------------------------------------------------------------------
bool NppD::Solve_Linear(double *Mat, double *rhs, bool resolve, int m, int nrhs)
{
  int info=0, l=m;
  char *cc="N";
  int *ipiv;
  ipiv = new int[l];

  if (!resolve) {
    (void) dgetrf_(&m, &m, Mat, &l, ipiv, &info);
    if (info < 0) cout << "ERROR dgetrf "<<info<<endl;
  }

  (void) dgetrs_(cc, &m, &nrhs, Mat, &l, ipiv, rhs, &m, &info);
  if (info < 0) cout << "ERROR dgetrs "<<info<<endl;
  
  
  //delete cc;
  delete [] ipiv;
}


//-----------------------------------------------------------------
// Function      : NppD::dphi_dlambda
// Purpose       : Use finite differences to calculate f(phi).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/05/05
//------------------------------------------------------------------
bool NppD::dphi_dlambda(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f)
{
  double delta = 1.0e-5; /*Finite Difference constant*/
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > phi;
  phi = createMember(f->space());

  App_Integrator->Integrate(xfinal,f,Tfinal,lambdafinal+delta);
  App_Integrator->Integrate(xfinal,phi,Tfinal,lambdafinal);
  Thyra::Vp_StV(&*f,-1.0,*phi);
  Thyra::Vt_S(&*f,1.0/delta);

  return true;
}
//-----------------------------------------------------------------
// Function      : NppD::Calculatenup
// Purpose       : Calculate the nu_p vector
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------

bool NppD::Calculatenup(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& nuq1,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& nuq2,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& nup1,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dphidlambda,
		       const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xi1,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& xi2,
		       double& deltalambda)
{
  double eps = 10e-5;

  // Dimension of the linear system depends on whether we are looking 
  // for a periodic solution.
  int LS_size = Unstable_Basis_Size+1;
    
  // First need a square MultiVector...
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > TempMV;
  TempMV = createMembers(xcurrent->space(),LS_size);
  double *LHS;
  double *RHS;
  LHS = new double[(LS_size)*(LS_size)];
  RHS = new double[LS_size];

  /*Upper left corner of the lhs */
  for (int j=0;j<Unstable_Basis_Size;j++)
    for (int i=0;i<Unstable_Basis_Size;i++)
      {
	LHS[i+j*(LS_size)]=Thyra::get_ele(*Re->col(j+1),i+1);
	if (i==j)
	  LHS[i+j*(LS_size)]+=-1.0;
      }

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp;
  temp = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > rightcol;
  rightcol = createMember(TempMV->domain());


  //Rightmost column of LHS.
  temp = MatVec(nuq2);
  Thyra::Vp_StV(&*temp,1.0,*xi2);
  Thyra::apply(*Vp,Thyra::TRANS,*temp,&*rightcol);

  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      LHS[i+(Unstable_Basis_Size)*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
    }
  
  //add a row.
  Thyra::apply(*Vp,Thyra::TRANS,*dvector,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      LHS[Unstable_Basis_Size+i*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
    }
  
  // and the corner
  LHS[(LS_size)*(LS_size)-1]=Thyra::dot(*dvector,*nuq2);
  
  // Now the Right Hand Side Vector
  temp = MatVec(nuq1);
  Thyra::Vp_StV(&*temp,1.0,*xi1);
  Thyra::apply(*Vp,Thyra::TRANS,*temp,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i]=-Thyra::get_ele(*rightcol,i+1);
    }

  Thyra::assign(&*temp,*nufinal);
  Thyra::Vp_StV(&*temp,1.0,*nuq1);
  RHS[Unstable_Basis_Size]=Thyra::dot(*dvector,*temp);

  Solve_Linear(LHS,RHS, false ,LS_size,1);
  for (int j=0;j<Unstable_Basis_Size;j++)
    {
      Thyra::set_ele(j+1,RHS[j],&*nup1);
    }

  deltalambda = RHS[Unstable_Basis_Size];

  delete [] LHS;
  delete [] RHS;

  return true;
}

//-----------------------------------------------------------------
// Function      : NppD::MT
// Purpose       : Calculate the M_T*u product
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
bool NppD::MT(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& MPartialVec,
	      const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& u,
	      const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v)
{
  // See Equation 35 in K.Engelborhgs Et. Al.

  double epsilon = 10e-8;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp;
  temp = createMember(u->space());

  temp = MatVec(u);
  Thyra::Vt_S(&*temp,epsilon);
  Thyra::Vp_StV(&*temp,1.0,*v);
  feval(MPartialVec,temp);
  feval(temp,v);
  Thyra::Vp_StV(&*MPartialVec,-1.0,*temp);
  Thyra::Vt_S(&*MPartialVec,1.0/epsilon);

  return true;
}
//-----------------------------------------------------------------
// Function      : NppD::MX
// Purpose       : Calculate the M_X*u product
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
bool NppD::MX(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& MPartialVec,
	      const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& u,
	      const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dx)
{
  // See Equation 35 in K.Engelborhgs Et. Al.
  double eps1 = 10e-8;
  double eps2 = 10e-8;

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp1;
  temp1 = createMember(u->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp2;
  temp2 = createMember(u->space());


  Thyra::assign(&*temp1,*xfinal);
  Thyra::Vp_StV(&*temp1,eps1,*u);
  Thyra::Vp_StV(&*temp1,eps2,*dx);
  
  App_Integrator->Integrate(temp1,MPartialVec,Tcurrent,lambdacurrent);

  Thyra::assign(&*temp1,*xfinal);
  Thyra::Vp_StV(&*temp1,eps2,*dx);
  App_Integrator->Integrate(temp1,temp2,Tcurrent,lambdacurrent);
  
  Thyra::Vp_StV(&*MPartialVec,-1.0,*temp2);

  Thyra::assign(&*temp1,*xfinal);
  Thyra::Vp_StV(&*temp1,eps1,*u);
  App_Integrator->Integrate(temp1,temp2,Tcurrent,lambdacurrent);
  
  Thyra::Vp_StV(&*MPartialVec,-1.0,*temp2);

  App_Integrator->Integrate(xfinal,temp2,Tcurrent,lambdacurrent);
  
  Thyra::Vp_StV(&*MPartialVec,1.0,*temp2);

  Thyra::Vt_S(&*MPartialVec,1.0/(eps1*eps2));
  return true;
}
//-----------------------------------------------------------------
// Function      : NppD::ML
// Purpose       : Calculate the M_L*u product
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/29/05
//------------------------------------------------------------------
bool NppD::ML(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& MPartialVec,
	      const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& u)
{
  // See Equation 35 in K.Engelborhgs Et. Al.
  double eps1 = 10e-8;
  double eps2 = 10e-8;

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp1;
  temp1 = createMember(u->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > temp2;
  temp2 = createMember(u->space());


  Thyra::assign(&*temp1,*xfinal);
  Thyra::Vp_StV(&*temp1,eps1,*u);
  
  App_Integrator->Integrate(temp1,MPartialVec,Tcurrent,lambdacurrent+eps2);

  //Thyra::assign(&*temp1,*xfinal);
  //Thyra::Vp_StV(&*temp1,eps2,*dx);
  App_Integrator->Integrate(xfinal,temp2,Tcurrent,lambdacurrent+eps2);
  
  Thyra::Vp_StV(&*MPartialVec,-1.0,*temp2);

  App_Integrator->Integrate(temp1,temp2,Tcurrent,lambdacurrent);
  Thyra::Vp_StV(&*MPartialVec,-1.0,*temp2);

  App_Integrator->Integrate(xfinal,temp2,Tcurrent,lambdacurrent);
  Thyra::Vp_StV(&*MPartialVec,1.0,*temp2);

  Thyra::Vt_S(&*MPartialVec,1.0/(eps1*eps2));
  return true;
}
