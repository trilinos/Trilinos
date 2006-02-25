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
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Npgs.hpp"

using namespace CAPO;

int Select(double *x, double *y);

//-----------------------------------------------------------------
// Function      : Npgs::Npgs
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Npgs::Npgs(const Teuchos::RefCountPtr<Parameter_List>& ParamList, 
	   const Teuchos::RefCountPtr<Integrator>& App_Int, 
	   const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x0, 
	   double lambda0, double T0) 
{
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
// Function      : Npgs::Initialize
// Purpose       : Create the initial subspace.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Initialize()
{
  unsigned int seedtime = time(NULL);
  Teuchos::ScalarTraits<Scalar>::seedrandom(seedtime);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > TempVector;
  TempVector = createMember(xcurrent->space());

  /*
    J. Simonis Error Report
    I should be able to seed the random number generator using the
    command Thyra::seed_randomize(seedtime); but this fails to compile
    so instead I must use the above seedrandom function.
  */

  for (int i=1;i<30+1;i++)
    {
      /*
      App_Integrator->Integrate(xcurrent,TempVector,(rand()%10)+20*Tfinal,lambdafinal);
      Thyra::Vp_S(&*TempVector,.01);
      Thyra::assign(&*Ve->col(i),*TempVector);
      */

      
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
      SchurDecomp(Se,Re);
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
// Function      : Npgs::Orthonormalize
// Purpose       : Orthonormalize the Vectors of a MultiVectorBase.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/23/05
//------------------------------------------------------------------
void Npgs::Orthonormalize(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Basis, int Number_of_Columns)
{
  double temp = 0;
  temp = sqrt( Thyra::dot( *(Basis->col(1)),*(Basis->col(1))));
  if (temp > 1.0e-12)
    Thyra::Vt_S( &*(Basis->col(1)), 1.0/temp );
  else
    cout <<"WARNING Npgs::Orthonormalize: norm=0; No scaling done" << endl;


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
    cout << "\n NPGS finished " << endl;
    //cout <<"\n Writing columns of Ve:" << endl;
  //Print(Ve);
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
// Creation Date : 06/20/05
//------------------------------------------------------------------
void Npgs::Predictor()
{

  // Eventually put into parameter list...
  double lambdamax = SolveParameters->get_lambda_max();
  double lambdamin = SolveParameters->get_lambda_min();

  Tstep = Tfinal-Tprevious;
  lambdastep = lambdafinal-lambdaprevious;
  Thyra::assign(&*xstep,*xfinal);
  Thyra::Vp_StV(&*xstep,-1.0,*xprevious);

  Tprevious = Tfinal;
  lambdaprevious = lambdafinal;
  Thyra::assign(&*xprevious,*xfinal);
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
    }
  else
    {
      Tcurrent = Tprevious + Tstep;
      lambdacurrent = lambdaprevious + lambdastep;

      Thyra::assign(&*xcurrent,*xprevious);
      Thyra::Vp_StV(&*xcurrent,1.0,*xstep);
    }

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
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > Npgs::MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y)
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
// Function      : Npgs::MatVecs
// Purpose       : returns M*Y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Npgs::MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Y)
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
  lambdastep = .5*lambdastep;
  Tstep = .5*Tstep;
  Thyra::Vt_S(&*xstep,.5);

  lambdacurrent = lambdaprevious + lambdastep;
  Tcurrent = Tprevious + Tstep;
  Thyra::assign(&*xcurrent,*xprevious);
  Thyra::Vp_StV(&*xcurrent,1.0,*xstep);

  cout << "Cutting the lambdastep in half." << endl;
  //abort();
}


//-----------------------------------------------------------------
// Function      : Npgs::InnerIteration
// Purpose       : The Npgs algorithm  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool Npgs::InnerIteration()
{
  bool converged=false;

  double eps = 10e-8;
  double *deltaT;
  deltaT = new double;
  double *deltalambda;
  deltalambda = new double;

  int Subspace_Size = 0;


  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > v;
  v = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > r;
  r = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dq;
  dq = createMember(xcurrent->space());

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
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dp;
	  dp = createMember(Vp->domain());

	  if ( SolveParameters->Arc_Length() )
	    {
	      ShermanMorrison(Vp,dq,dp,Re,v,finit,r,*deltaT,*deltalambda);
	      Thyra::Vp_StV(&*xfinal,1.0,*dq);
	      Thyra::apply(*Vp,Thyra::NOTRANS,*dp,&*TempVector);
	      Thyra::Vp_StV(&*xfinal,1.0,*TempVector);
	      Tfinal+=*deltaT;
	      lambdafinal+=*deltalambda;

	    }
	  else
	    {	  
	      Calculatedq(Vp,dq,r);
	      Calculatedp(Vp,dq,dp,Re,v,finit,r,*deltaT);
	      Thyra::Vp_StV(&*xfinal,1.0,*dq);
	      Thyra::apply(*Vp,Thyra::NOTRANS,*dp,&*TempVector);
	      Thyra::Vp_StV(&*xfinal,1.0,*TempVector);
	      Tfinal+=*deltaT;
	    }
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
	      dq = MatVec(TempVector);
	      Thyra::Vp_StV(&*dq,1.0,*r);
	      Thyra::Vp_StV(&*xfinal,1.0,*dq);

	    }
	  else // Not a periodic problem.
	    {
	      if ( SolveParameters->Arc_Length() )
		{
		  TempVector = MatVec(r);
		  Thyra::Vp_StV(&*TempVector,1.0,*r);
		  dq = MatVec(TempVector);
		  Thyra::Vp_StV(&*dq,1.0,*r);

		  //essentially set deltaT = 0 by not changing Tfinal.
		  dphi_dlambda(TempVector);
		  *deltalambda = ( Thyra::dot(*TempVector,*dq) )/lambdastep;
		  lambdafinal+=*deltalambda;
		  Thyra::Vp_StV(&*xfinal,1.0,*dq);
		}
	      else
		{
		  // Just set deltaT = 0, or really leave Tfinal alone.
		  // Still update x(0), here assume l=2 (still hardcoded...)
		  TempVector = MatVec(r);
		  Thyra::Vp_StV(&*TempVector,1.0,*r);
		  dq = MatVec(TempVector);
		  Thyra::Vp_StV(&*dq,1.0,*r);
		  Thyra::Vp_StV(&*xfinal,1.0,*dq);
		}
	    }
	  /*
	  dq = MatVec(r);
	  Thyra::Vp_StV(&*dq,1.0,*r);

	  App_Integrator->Integrate(xfinal,v,eps,lambdafinal);
	  Thyra::Vp_StV(&*v,-1.0,*xfinal);
	  Thyra::Vt_S(&*v,1.0/eps);
	  double lhs = Thyra::dot(*finit,*v);

	  Thyra::assign(&*TempVector,*xfinal);
	  Thyra::Vp_StV(&*TempVector,-1.0,*xcurrent);
	  Thyra::Vp_StV(&*TempVector,1.0,*dq);
	  double rhs = Thyra::dot(*finit,*TempVector);

	  *deltaT = -rhs/lhs;

	  Tfinal +=(*deltaT);
	  Thyra::Vp_StV(&*xfinal,1.0,*dq);
	  */
	}

      // If my subspace has shrunk at any step, I need to replace the vectors
      // in Ve with new randomized vectors.
      /*
      for (int k = Unstable_Basis_Size+SolveParameters->get_NumberXtraVecsSubspace()+1;k<31;k++)
	{
	  
	  App_Integrator->Integrate(xcurrent,TempVector,(rand()%10)+20*Tfinal,lambdafinal);
	  Thyra::Vp_S(&*TempVector,.01);
	  Thyra::assign(&*Ve->col(k),*TempVector);
	  
	  
	  Thyra::randomize(0.1,1.1,&*Ve->col(k));
	  App_Integrator->Integrate(Ve->col(k),TempVector,5*Tfinal,lambdafinal);
	  Thyra::assign(&*Ve->col(k),*TempVector);
	  
	}
      */
      Orthonormalize(Ve,Unstable_Basis_Size+SolveParameters->get_NumberXtraVecsSubspace());
      App_Integrator->Integrate(xfinal,v,Tfinal,lambdafinal);
      Thyra::assign(&*r,*v); 
      Thyra::Vp_StV(&*r,-1.0,*xfinal);

      converged = Converged(xfinal,v);


      iter++;
    }//end while

  delete deltaT;
  delete deltalambda;

  return converged;

}


//-----------------------------------------------------------------
// Function      : Npgs::SchurDecomposition
// Purpose       : One Multivector comes in, Schur Decomposition
//                 is performed on it, and the vectors are returned
//                 along with a multivector containing the scaling factors.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
void Npgs::SchurDecomp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re)
{
  /* This function performs a Schur Decomposition.
     On return Mat contains the Upper Triangular matrix with 
     eigenvalues on the diagonal, and V contains the basis
     vectors.  This function also requires the existence of
     function SELECT which is called in when ordering the 
     eigenvalues along the diagonal.
  */

  Teuchos::LAPACK<int,double> Tlapack;

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
  
  pt2Select = &Select;

  Tlapack.GEES(*cn, *cs, pt2Select, m, se, LDA, &sdim, 
	 wr, wi, re, LDA, work, lwork, bwork, &info);

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
// Function      : Npgs::Print
// Purpose       : Print a MultiVector
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Npgs::Print(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Printme)
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
// Function      : Npgs::Select
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
     an int rather than a bool.  Thus, true=0, false=1.  For now
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
// Function      : Npgs::Converged
// Purpose       : Checks to see if the norm(x-y)<tol
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool Npgs::Converged(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,
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
// Function      : Npgs::SubspaceIterations
// Purpose       : See Lust et.al.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void Npgs::SubspaceIterations(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We, const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re)
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
      SchurDecomp(Se,Re);
      if (i<(SolveParameters->get_SubspaceIterations()-1))
	{
	  Thyra::apply(*We,Thyra::NOTRANS,*Se,&*Ve_pe);
	  Orthonormalize(Ve_pe,Subspace_Size);
	}
    }
  Orthonormalize(Ve_pe,Subspace_Size);
}

//-----------------------------------------------------------------
// Function      : Npgs::Calculatedq
// Purpose       : Perform the Picard iterations
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void Npgs::Calculatedq(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp
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
// Function      : Npgs::ComputeVp
// Purpose       : Extract Vp from Ve and multiply it by Se.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Npgs::ComputeVp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,
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
// Function      : Npgs::UpdateVe
// Purpose       : Update the Basis
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Npgs::UpdateVe(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We,
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
// Function      : Npgs::Calcuatedp
// Purpose       : Find the dp vector and deltaT
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Npgs::Calculatedp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dp,
		       const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
		       double& deltaT)
{
  double eps = 10e-5;


  // Dimension of the linear system depends on whether we are looking 
  // for a periodic solution.
  int LS_size = Unstable_Basis_Size;

  if (SolveParameters->Periodic())
    LS_size++;
    
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
  
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > fphi;
  fphi = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > rightcol;
  rightcol = createMember(TempMV->domain());
  
  //Rightmost column of LHS.
  if (SolveParameters->Periodic())
    {
      // If periodic solution we add a column and a row...
      dphi_dt(fphi);

      Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol);
 
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
      LHS[(Unstable_Basis_Size+1)*(Unstable_Basis_Size+1)-1]=0.0;

    }

  //Have Built the Left hand side, now need the right hand side.

  fphi = MatVec(dq);
  Thyra::Vp_StV(&*fphi,1.0,*r);
  // I use the rightcol vector for storage, just so I don't have to
  // create an extra temporary vector...
  Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol);
  // rightcol now contains the first p elements 

  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i]=-Thyra::get_ele(*rightcol,i+1);
    }

  if (SolveParameters->Periodic())
    {
      Thyra::assign(&*fphi,*xfinal);
      Thyra::Vp_StV(&*fphi,-1.0,*xinit);
      Thyra::Vp_StV(&*fphi,1.0,*dq);
      
      RHS[Unstable_Basis_Size]=-Thyra::dot(*finit,*fphi);
    }

  Solve_Linear(LHS,RHS, false ,LS_size,1);
  for (int j=0;j<Unstable_Basis_Size;j++)
    {
      Thyra::set_ele(j+1,RHS[j],&*dp);
    }

  if (SolveParameters->Periodic())
    deltaT = RHS[Unstable_Basis_Size];

  delete [] LHS;
  delete [] RHS;

  return true;
}


//-----------------------------------------------------------------
// Function      : Npgs::dphi_dt
// Purpose       : Use finite differences to calculate f(phi).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Npgs::dphi_dt(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f)
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
// Function      : Npgs::Solve_Linear
// Purpose       : Call LAPACK routine to solve a linear system
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/20/05
//------------------------------------------------------------------
bool Npgs::Solve_Linear(double *Mat, double *rhs, bool resolve, int m, int nrhs)
{
  Teuchos::LAPACK<int,double> Tlapack;

  int info=0, l=m;
  char *cc="N";
  int *ipiv;
  ipiv = new int[l];

  if (!resolve) {
    Tlapack.GETRF(m, m, Mat, l, ipiv, &info);
    if (info < 0) cout << "ERROR dgetrf "<<info<<endl;
  }

  Tlapack.GETRS(*cc, m, nrhs, Mat, l, ipiv, rhs, m, &info);
  if (info < 0) cout << "ERROR dgetrs "<<info<<endl;
  
  
  //delete cc;
  delete [] ipiv;
}


//-----------------------------------------------------------------
// Function      : Npgs::ShermanMorrison
// Purpose       : Find the dp vector and deltaT
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/05/05
//------------------------------------------------------------------
bool Npgs::ShermanMorrison(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dp,
		       const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
		       const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
		       double& deltaT, double& deltalambda)
{
  double eps = 10e-5;

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > fphi;
  fphi = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > templong;
  templong = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dq1;
  dq1 = createMember(xcurrent->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dq2;
  dq2 = createMember(xcurrent->space());

  // Need to caluclate the two different dq's for the two different
  // linear systems.
  dphi_dlambda(fphi);
  Thyra::Vt_S(&*fphi,-1.0);
  Calculatedq(Vp,dq1,fphi);
  Calculatedq(Vp,dq2,r);

  // Dimension of the linear system depends on whether we are looking 
  // for a periodic solution and/or performing psuedo-arclength
  // continuation.
  int LS_size = Unstable_Basis_Size+1;

  if (SolveParameters->Periodic())
    LS_size++;
    
  // First need a square MultiVector...
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > TempMV;
  TempMV = createMembers(xcurrent->space(),LS_size);
  double *LHS;
  double *RHS;
  LHS = new double[(LS_size)*(LS_size)];
  RHS = new double[2*LS_size];

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > rightcol;
  rightcol = createMember(TempMV->domain());

  /*Upper left corner of the lhs */
  for (int j=0;j<Unstable_Basis_Size;j++)
    for (int i=0;i<Unstable_Basis_Size;i++)
      {
	LHS[i+j*(LS_size)]=Thyra::get_ele(*Re->col(j+1),i+1);
	if (i==j)
	  LHS[i+j*(LS_size)]+=-1.0;
      }

  /* and the first p elements of RHS */
  fphi = MatVec(dq1);
  Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol);
  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i]=-Thyra::get_ele(*rightcol,i+1);
    }
  Thyra::assign(&*templong,*dq2);
  fphi = MatVec(templong);
  Thyra::Vp_StV(&*fphi,1.0,*r);


  Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol);

  for (int i=0;i<Unstable_Basis_Size;i++)
    {
      RHS[i+LS_size]=-Thyra::get_ele(*rightcol,i+1);
    }

  if ( SolveParameters->Periodic() )
    {
      // Add the final two rows and columns onto LHS

      // Add column 1
      dphi_dt(fphi);
      Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol); 
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[i+(Unstable_Basis_Size)*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}

      //add row 1.

      Thyra::apply(*Vp,Thyra::TRANS,*finit,&*rightcol);
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[Unstable_Basis_Size+i*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}

      // and corner 1
      LHS[(Unstable_Basis_Size+1)*(LS_size)-2]=0.0;

      // Add column 2
      dphi_dlambda(fphi);
      Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol); 
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[i+(Unstable_Basis_Size+1)*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}

      //add row 2.
      // Assume psuedo arclength condition given by eqn. 4.6 in Lust et. al.
      Thyra::apply(*Vp,Thyra::TRANS,*xstep,&*rightcol);
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[Unstable_Basis_Size+1+i*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}
      Thyra::assign(&*fphi,*xfinal);
      Thyra::Vp_StV(&*fphi,-1.0,*xcurrent);

      /*
      RHS[2*LS_size-1]=-(Thyra::dot(*xstep,*fphi)+(Tfinal-Tcurrent)*(Tstep)+(lambdastep)*(lambdafinal-lambdacurrent))-Thyra::dot(*xstep,*dq2);
      */
      RHS[2*LS_size-1]=-Thyra::dot(*xstep,*dq2);
      RHS[LS_size-1]=-Thyra::dot(*xstep,*dq1);

      // and corners 2 and 3
      LHS[(LS_size)*(Unstable_Basis_Size+1)-1]=Tstep;
      LHS[(LS_size)*(LS_size)-1]=(lambdastep);
      LHS[(LS_size)*(LS_size)-2]=0.0;

      Thyra::assign(&*fphi,*xfinal);
      Thyra::Vp_StV(&*fphi,-1.0,*xinit);
      Thyra::Vp_StV(&*fphi,1.0,*dq2);
      RHS[2*LS_size-2]=-Thyra::dot(*finit,*fphi);

      RHS[LS_size-2]=-Thyra::dot(*finit,*dq1);

    }
  else
    {
      // Add column 1
      dphi_dlambda(fphi);
      Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol); 
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[i+(Unstable_Basis_Size)*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}

      //add row 1.
      // Assume psuedo arclength condition given by eqn. 4.6 in Lust et. al.
      Thyra::assign(&*fphi,*xcurrent);
      Thyra::Vp_StV(&*fphi,-1.0,*xprevious);
      Thyra::apply(*Vp,Thyra::TRANS,*fphi,&*rightcol);
      for (int i=0;i<Unstable_Basis_Size;i++)
	{
	  LHS[Unstable_Basis_Size+i*(LS_size)]=Thyra::get_ele(*rightcol,i+1);
	}
      RHS[LS_size-1]=-Thyra::dot(*fphi,*dq1);

      // and corner
      LHS[(LS_size)*(LS_size)-1]=lambdastep;

      // Last element of the right hand sides...
      /*
      Thyra::assign(&*fphi,*xfinal);
      Thyra::Vp_StV(&*fphi,-1.0,*xcurrent);
      RHS[2*LS_size-1]=-Thyra::dot(*fphi,*xstep)-(Tfinal-Tcurrent)*Tstep-
	(lambdafinal-lambdacurrent)*lambdastep-Thyra::dot(*xstep,*dq);
      */
      RHS[2*LS_size-1]=-Thyra::dot(*fphi,*dq2);
    }


  Solve_Linear(LHS,RHS, false ,LS_size,2);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > omega;
  omega = createMember(dp->space());
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > nu;
  nu = createMember(dp->space());

  deltalambda = (RHS[2*LS_size-1])/(1.0+RHS[LS_size-1]);

  for (int j=0;j<Unstable_Basis_Size;j++)
    {
      Thyra::set_ele(j+1,RHS[LS_size+j]-(deltalambda)*RHS[j],&*dp);
    }

  if (SolveParameters->Periodic())
    {
      deltaT = RHS[2*LS_size-2]-(deltalambda)*RHS[LS_size-2];
    }
  Thyra::assign(&*dq,*dq2);
  Thyra::Vp_StV(&*dq,-deltalambda,*dq1);
  delete [] LHS;
  delete [] RHS;

  return true;
}


//-----------------------------------------------------------------
// Function      : Npgs::dphi_dlambda
// Purpose       : Use finite differences to calculate f(phi).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/05/05
//------------------------------------------------------------------
bool Npgs::dphi_dlambda(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f)
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

