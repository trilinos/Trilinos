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
File:      Capo_Rpm.cpp
Purpose:   The Newton--Picard Gauss-Seidel Solver for Capo.
Date:      6-30-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Rpm.hpp"

using namespace CAPO;

//-----------------------------------------------------------------
// Function      : Rpm::Rpm
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
Rpm::Rpm(const Teuchos::RefCountPtr<Parameter_List>& ParamList, \
	 const Teuchos::RefCountPtr<Integrator>& App_Int, \
	 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x0,\
	 double lambda0, double T0)
{
  xcurrent = createMember(x0->space());
  xfinal = createMember(x0->space());
  Thyra::assign(&*xcurrent, *x0); 
  Thyra::assign(&*xfinal, *x0);

  lambdacurrent = lambda0;
  lambdafinal = lambda0;

  Tcurrent = T0;
  Tfinal = T0;

  iter = 0;
  Unstable_Basis_Size = 0;

  App_Integrator = App_Int;
  SolveParameters = ParamList;

  // Use the VectorBase xcurrent to create a MultiVector Base.
  Ve = Thyra::createMembers(xcurrent->space(),30);

}

//-----------------------------------------------------------------
// Function      : Rpm::Initialize
// Purpose       : Create the initial subspace.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Rpm::Initialize()
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
// Function      : Rpm::Orthonormalize
// Purpose       : Orthonormalize the Vectors of a MultiVectorBase.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/23/05
//------------------------------------------------------------------
void Rpm::Orthonormalize(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Basis, int Number_of_Columns)
{
  double temp = 0;
  temp = sqrt( Thyra::dot( *(Basis->col(1)),*(Basis->col(1))));
  if (temp > 1.0e-12)
    Thyra::Vt_S( &*(Basis->col(1)), 1.0/temp );
  else
    cout <<"WARNING Rpm::Orthonormalize: norm=0; No scaling done" << endl;


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
	cout <<"WARNING Rpm::Orthonormalize: norm=0; No scaling done" << endl;
    }
}


//-----------------------------------------------------------------
// Function      : Rpm::Finish
// Purpose       : Print the Basis Matrix
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Rpm::Finish()
{
  if ((*SolveParameters).get_printproc() >0)
    cout << "\n NPGS finished " << endl;
    //cout <<"\n Writing columns of Ve:" << endl;
  //Print(Ve);
}
//-----------------------------------------------------------------
// Function      : Rpm::InnerFunctions
// Purpose       : Rpm does not require anything to be done 
//                 between inner iterations.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Rpm::InnerFunctions()
{
}
//-----------------------------------------------------------------
// Function      : Rpm::Predictor
// Purpose       : Currently increments lambda, and sets the new
//                 solution guesses to the current solution.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/20/05
//------------------------------------------------------------------
void Rpm::Predictor(double& StepSize,double& PrevStepSize)
{
  lambdacurrent = lambdafinal+(*SolveParameters).get_lambda_stepsize();
  Tcurrent = Tfinal;
  Thyra::assign(&*xcurrent,*xfinal);
}
//-----------------------------------------------------------------
// Function      : Rpm::Get_Tfinal
// Purpose       : Returns Tfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& Rpm::Get_Tfinal()
{
  return Tfinal;
}


//-----------------------------------------------------------------
// Function      : Rpm::Get_lambdafinal
// Purpose       : Returns lambdafinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
double& Rpm::Get_lambdafinal()
{
  return lambdafinal;
}


//-----------------------------------------------------------------
// Function      : Rpm::Get_xfinal
// Purpose       : Returns xfinal
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Rpm::Get_xfinal()
{
  return xfinal;
}


//-----------------------------------------------------------------
// Function      : Rpm::MatVec
// Purpose       : returns M*y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > Rpm::MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y)
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
// Function      : Rpm::MatVecs
// Purpose       : returns M*Y, a matrix vector product obtained
//                 using finite differences.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Rpm::MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Y)
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
// Function      : Rpm::IterationFailure
// Purpose       : If the inner iteration function fails to find a
//                 solution, then the code must exit.  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Rpm::IterationFailure()
{
  cout << "Rpm::IterationFailure Failed to find a solution." << endl;
  cout << "I should find a better way of doing this..." << endl;
  abort();
}


//-----------------------------------------------------------------
// Function      : Rpm::InnerIteration
// Purpose       : The Rpm algorithm  
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool Rpm::InnerIteration()
{
  bool converged=false;

  double eps = 10e-8;
  double *deltaT;
  deltaT = new double;

  int Subspace_Size = 0;

  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > finit;
  finit = createMember(xcurrent->space());

  App_Integrator->Integrate(xcurrent,finit,eps,lambdacurrent);
  Thyra::Vp_StV(&*finit,-1.0,*xcurrent);
  Thyra::Vt_S(&*finit,1.0/eps);

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
      //Orthonormalize(Ve,Subspace_Size);
      SubspaceIterations(Se,We,Re);

      if (Unstable_Basis_Size > 0) // Need to include a Newton-Step.
	{
	  // Declarations (these change size with every go-around)
	  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Vp;
	  Vp = createMembers(xcurrent->space(),Unstable_Basis_Size);

	  ComputeVp(Se,Vp);
	  Calculatedq(Vp,dq,r);
	  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > dp;
	  dp = createMember(Vp->domain());
	  Calculatedp(Vp,dq,dp,Re,v,finit,r,*deltaT);
	  Thyra::Vp_StV(&*xfinal,1.0,*dq);
	  Thyra::apply(*Vp,Thyra::NOTRANS,*dp,&*TempVector);
	  Thyra::Vp_StV(&*xfinal,1.0,*TempVector);
	  Tfinal+=*deltaT;
	  UpdateVe(We,Se);
	}
      else
	{
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
	}


      Orthonormalize(Ve,Unstable_Basis_Size+SolveParameters->get_NumberXtraVecsSubspace());
      App_Integrator->Integrate(xfinal,v,Tfinal,lambdafinal);
      Thyra::assign(&*r,*v); 
      Thyra::Vp_StV(&*r,-1.0,*xfinal);

      converged = Converged(xfinal,v);


      iter++;
    }//end while

  return true;

}


//-----------------------------------------------------------------
// Function      : Rpm::SchurDecomposition
// Purpose       : One Multivector comes in, Schur Decomposition
//                 is performed on it, and the vectors are returned
//                 along with a multivector containing the scaling factors.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
void Rpm::SchurDecomp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re)
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
  
  pt2Select = &Select;

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
// Function      : Rpm::Print
// Purpose       : Print a MultiVector
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/13/05
//------------------------------------------------------------------
void Rpm::Print(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Printme)
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
// Function      : Rpm::Select
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
// Function      : Rpm::Converged
// Purpose       : Checks to see if the norm(x-y)<tol
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/15/05
//------------------------------------------------------------------
bool Rpm::Converged(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,
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
// Function      : Rpm::SubspaceIterations
// Purpose       : See Lust et.al.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void Rpm::SubspaceIterations(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We, const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re)
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
// Function      : Rpm::Calculatedq
// Purpose       : Perform the Picard iterations
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/16/05
//------------------------------------------------------------------
void Rpm::Calculatedq(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp
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
}


//-----------------------------------------------------------------
// Function      : Rpm::ComputeVp
// Purpose       : Extract Vp from Ve and multiply it by Se.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Rpm::ComputeVp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,
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
// Function      : Rpm::UpdateVe
// Purpose       : Update the Basis
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Rpm::UpdateVe(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We,
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
// Function      : Rpm::Calcuatedp
// Purpose       : Find the dp vector and deltaT
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Rpm::Calculatedp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
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
  // for a periodic solution and/or performing psuedo-arclength
  // continuation.
  int LS_size = Unstable_Basis_Size;

  if (SolveParameters->Periodic())
    LS_size++;
  if (SolveParameters->Arc_Length())
    LS_size++;


  // First need a p+1 by p+1 MultiVector...
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > TempMV;
  TempMV = createMembers(xcurrent->space(),LS_size);
  double *LHS;
  double *RHS;
  LHS = new double[(LS_size)*(LS_Size)];
  RHS = new double[LS_size];

  /*Upper left corner of the lhs */
  for (int j=0;j<Unstable_Basis_Size;j++)
    for (int i=0;i<Unstable_Basis_Size;i++)
      {
	LHS[i+j*(LS_size)]=Thyra::get_ele(*Re->col(j+1),i+1);
	if (i==j)
	  LHS[i+j*(LS_size)]+=-1.0;
      }


  //Rightmost column of LHS.
  if (SolveParameters->Periodic())
    {
      // If periodic solution we add a column and a row...
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > fphi;
      fphi = createMember(xcurrent->space());
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > rightcol;
      rightcol = createMember(TempMV->domain());
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
      
      App_Integrator->Integrate(xfinal,fphi,eps,lambdacurrent);
      Thyra::Vp_StV(&*fphi,-1.0,*xfinal);
      Thyra::Vt_S(&*fphi,1.0/eps);
      
      LHS[(Unstable_Basis_Size+1)*(Unstable_Basis_Size+1)-1]=Thyra::dot(*finit,*fphi);

    }
  // When psuedo-arclength continuation is added in, I'll need another row and 
  // column here, but for now I'll ignore these.

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
      Thyra::Vp_StV(&*fphi,-1.0,*xcurrent);
      Thyra::Vp_StV(&*fphi,1.0,*dq);
      
      RHS[Unstable_Basis_Size]=-Thyra::dot(*finit,*fphi);
    }

  Solve_Linear(LHS,RHS, false ,LS_size);
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
// Function      : Rpm::dphi_dt
// Purpose       : Use finite differences to calculate f(phi).
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/17/05
//------------------------------------------------------------------
bool Rpm::dphi_dt(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f)
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
// Function      : Rpm::Solve_Linear
// Purpose       : Call LAPACK routine to solve a linear system
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/20/05
//------------------------------------------------------------------
bool Rpm::Solve_Linear(double *Mat, double *rhs, bool resolve, int m)
{
  int info=0, nrhs=1, l=m;
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





}
