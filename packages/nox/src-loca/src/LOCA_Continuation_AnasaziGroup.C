// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Continuation_AnasaziGroup.H"
#include "NOX_Parameter_List.H"
#include "LOCA_Utils.H"
#include "LOCA_TimeDependent_AbstractGroup.H"

#ifdef HAVE_LOCA_ANASAZI
#include "AnasaziLOCAInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#endif

LOCA::Continuation::AnasaziGroup::AnasaziGroup() :
  eigenvalCounter(0),
  hasMassMatrix(false),
  cayley(false),
  sigma(0.0),
  mu(0.0),
  shift_(0.0)
{
}

LOCA::Continuation::AnasaziGroup::AnasaziGroup(
			   const LOCA::Continuation::AnasaziGroup& source, 
			   NOX::CopyType type) :
  eigenvalCounter(source.eigenvalCounter),
  hasMassMatrix(source.hasMassMatrix),
  cayley(source.cayley),
  sigma(source.sigma),
  mu(source.mu),
  shift_(source.shift_)
{
}

LOCA::Continuation::AnasaziGroup::~AnasaziGroup()
{
}

LOCA::Continuation::AnasaziGroup&
LOCA::Continuation::AnasaziGroup::operator=(
			     const LOCA::Continuation::AnasaziGroup& source)
{
  if (this != &source) {
    eigenvalCounter = source.eigenvalCounter;
  }

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AnasaziGroup:: computeEigenvalues(
					       NOX::Parameter::List& params)
{
#ifdef HAVE_LOCA_ANASAZI

  NOX::Parameter::List& aList = params.sublist("LOCA").sublist("Stepper").sublist("Anasazi");
  // The following are Parameter for the Anasazi Eigensolver
  int blksz =   aList.getParameter("Block Size", 1);       //  The block size
  int length =  aList.getParameter("Arnoldi Size", 30);      //  The maximum length of the Arnoldi factorization
  int nev =     aList.getParameter("NEV", 4);         //  The number of requested eigenvalues
  double tol =  aList.getParameter("Tol", 1.0e-7); //  Tolerance for the converged eigenvalues
  int step =    aList.getParameter("Convergence Check", 1);        //  This checks convergence every so many steps
  int restart = aList.getParameter("Restarts",1);    //  This is the number of restarts allowed
  int freq =    aList.getParameter("Frequency",1);    // How often to recalculate eigenvalues
  int debug =   aList.getParameter("Debug Level",1);  // Anasazi Debug level
  string which= aList.getParameter("Sorting Order","LM");   //  Which eigenvalues are of interest.
  hasMassMatrix = aList.getParameter("Mass Matrix",false);  //Is there a mass matrix
  cayley = aList.getParameter("Cayley Transformation",false); // Do Cayley transformation
  sigma = aList.getParameter("Cayley Pole", 0.0);  // Cayley pole
  mu = aList.getParameter("Cayley Zero", 0.0); // Cayley zero
  shift_ = aList.getParameter("Shift",0.0);   // Shift for shift-and-invert

  // Check if eigenvalues are requested this continuation step
  if (eigenvalCounter++%freq != 0) {
    if (Utils::doPrint(Utils::StepperIteration)) 
      cout <<"\tAnasazi Eigensolver not requested this continuation step." 
	   << endl;
    return NOX::Abstract::Group::Ok;
  }
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\n" << Utils::fill(64,'=')
         << "\nAnasazi Eigensolver starting with block size " << blksz
	 << "\n" << endl;
  }

  // Create updated Jacobian matrix
  computeJacobian();

   // Get reference to solution vector to clone
  const NOX::Abstract::Vector& xVector = getX();

  // Create the operator and initial vector
  Anasazi::LOCAMat<double> Amat( params, *this );
  Anasazi::LOCAVec<double> ivec( xVector, blksz );
  ivec.MvRandom();

  // Apply J^-1*M to random vector before eigensolve begins

  /* 
      NOX::Abstract::Vector *noxivec =  0;
      NOX::Abstract::Vector *ivecclone = xVector.clone(NOX::ShapeCopy);
      NOX::Abstract::Vector *intermed = xVector.clone(NOX::ShapeCopy);

     for (int j=0; j<blksz; j++){
    noxivec = ivec.GetNOXVector( j );
    ivecclone->update(1.0,*noxivec);
    if(hasMassMatrix) {
    const LOCA::TimeDependent::AbstractGroup& timegroup = 
    dynamic_cast<const LOCA::TimeDependent::AbstractGroup&>(*this);
    timegroup.applyMassMatrix(*ivecclone,*intermed);
    applyJacobianInverse(params,*intermed,*noxivec);
    }
    else {
      applyJacobianInverse(params,*ivecclone,*noxivec);
    }
   }
*/

  // Create an instance of the eigenproblem
  Anasazi::Eigenproblem<double> LOCAProblem( &Amat, &ivec );

  // Initialize the solver
  Anasazi::BlockArnoldi<double> LOCABlockArnoldi(LOCAProblem, tol, nev, length,
						 blksz, which, step, restart );

  // Print out debugging information on single proc
  LOCABlockArnoldi.setDebugLevel(debug);

  // Solve the problem to the specified tolerance
  LOCABlockArnoldi.solve();

  // Look at the solutions once if debug=0
  if (Utils::doPrint(Utils::StepperIteration)) 
    if (debug == 0) LOCABlockArnoldi.currentStatus();

  // Obtain the eigenvalues / eigenvectors
  int narn =  length; 
  double * evalr = LOCABlockArnoldi.getEvals(narn);  // narn modified
  double * evali = LOCABlockArnoldi.getiEvals(narn); // to within [nev,length]

  if (Utils::doPrint(Utils::StepperIteration)) {
    if(cayley){
          cout<<"Untransformed eigenvalues (since the operator was Cayley transformation)"<<endl;
    }
    else{
      if(shift_ == 0.0){
       cout<<"Untransformed eigenvalues (since the operator was the Jacobian inverse)"<<endl;
      }
      else{
             cout<<"Untransformed eigenvalues (since the operator was shift-and-invert)"<<endl;
      }
    }
  }
  
  // Obtain the eigenvectors
  Anasazi::LOCAVec<double> evecR( xVector, nev );
  LOCABlockArnoldi.getEvecs( evecR );
  Anasazi::LOCAVec<double> evecI( xVector, nev );
  LOCABlockArnoldi.getiEvecs( evecI );

  // Create updated Jacobian matrix
  computeJacobian();

  // Create some temporary vectors
  NOX::Abstract::Vector *r_evec = 0;
  NOX::Abstract::Vector *i_evec = 0;
  NOX::Abstract::Vector *tempvecr = xVector.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *tempveci = xVector.clone(NOX::ShapeCopy);
  NOX::Abstract::Group::ReturnType res;
  double realpart, imagpart; 
  for (int i=0; i<nev; i++) {

    // Computes z^T(Jz) for each eigenvector z

    r_evec = evecR.GetNOXVector( i );
    i_evec = evecI.GetNOXVector( i );
    res = applyJacobian(*r_evec, *tempvecr);
    res = applyJacobian(*i_evec, *tempveci);
    realpart = r_evec->dot(*tempvecr)+i_evec->dot(*tempveci);
    imagpart = r_evec->dot(*tempveci)-i_evec->dot(*tempvecr);

    // If you have mass matrix, compute z^T(Mz) for each eigenvector
    // and then [z^T(Jz)]/[z^T(Mz)]

    if (hasMassMatrix) {
      double rquot, iquot, rtemp, itemp, magn;
      NOX::Abstract::Vector *r_mass = xVector.clone(NOX::ShapeCopy);
      NOX::Abstract::Vector *i_mass = xVector.clone(NOX::ShapeCopy);
      const LOCA::TimeDependent::AbstractGroup& timeGroup = 
      dynamic_cast<const LOCA::TimeDependent::AbstractGroup&>(*this);
      timeGroup.applyMassMatrix(*r_evec, *r_mass);
      timeGroup.applyMassMatrix(*i_evec, *i_mass);
      rquot = r_evec->dot(*r_mass)+i_evec->dot(*i_mass);
      iquot = r_evec->dot(*i_mass)-i_evec->dot(*r_mass);
      magn = rquot*rquot + iquot*iquot;
      rtemp = realpart;
      itemp = imagpart;
      realpart = (rtemp*rquot + itemp*iquot)/magn;
      imagpart = (itemp*rquot - rtemp*iquot)/magn;
      delete r_mass;
      delete i_mass;
    }

    // Print out eigenvalue and Rayleigh quotient residual

    if (Utils::doPrint(Utils::StepperIteration)) {
      if(cayley) {
	double omag = (evalr[i]-1)*(evalr[i]-1)+evali[i]*evali[i];           
        double mag=evalr[i]*evalr[i]+evali[i]*evali[i];             
        double eigreal,eigimag;

        eigreal = (sigma*mag - (sigma + mu)*evalr[i] + mu)/omag;
        eigimag = (mu - sigma)*evali[i]/omag;

      cout<<"Eigenvalue "<<i<<" : "<<LOCA::Utils::sci(eigreal)<<"  "
	  <<LOCA::Utils::sci(eigimag)<<" i    :  RQresid "
          << LOCA::Utils::sci(fabs(eigreal - realpart)) <<"  "
	  << LOCA::Utils::sci(fabs(eigimag - imagpart))<<" i"<<endl;
      }
      else {
      double mag=evalr[i]*evalr[i]+evali[i]*evali[i];
      cout<<"Eigenvalue "<<i<<" : "<<LOCA::Utils::sci(evalr[i]/mag+shift_)<<"  "
	  <<LOCA::Utils::sci(-evali[i]/mag)<<" i    :  RQresid "
          << LOCA::Utils::sci(fabs(evalr[i]/mag + shift_ - realpart)) <<"  "
	  << LOCA::Utils::sci(fabs(-evali[i]/mag - imagpart))<<" i"<<endl;
      }  
    }
  }  

  // Print out remaining eigenvalue approximations from nev to final arnoldi size
  if (Utils::doPrint(Utils::StepperIteration) && narn>nev) {
    cout << "~~~~~~~ remaining eigenvalue approximations ~~~~~~~~~~~~" << endl;
    for (int i=nev; i<narn; i++) {
      if(cayley) {
	double omag = (evalr[i]-1)*(evalr[i]-1)+evali[i]*evali[i];           
        double mag=evalr[i]*evalr[i]+evali[i]*evali[i];             
        double eigreal,eigimag;

        eigreal = (sigma*mag - (sigma + mu)*evalr[i] + mu)/omag;
        eigimag = (mu - sigma)*evali[i]/omag;

      cout<<"Eigenvalue "<<i<<" : "<<LOCA::Utils::sci(eigreal)<<"  "
	  <<LOCA::Utils::sci(eigimag)<<" i"<<endl;
      }
      else {
        double mag=evalr[i]*evalr[i]+evali[i]*evali[i];
        cout<<"Eigenvalue "<<i<<" : "<<LOCA::Utils::sci(evalr[i]/mag+shift_)<<"  "
	    <<LOCA::Utils::sci(-evali[i]/mag)<<" i"<<endl;
      }
    }
  }

  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\nAnasazi Eigensolver finished.\n" 
         << Utils::fill(64,'=') << "\n" << endl;
  }

  delete [] evalr;
  delete [] evali;
  delete tempvecr;
  delete tempveci;

  return NOX::Abstract::Group::Ok;
#else
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\nWarning: LOCA::Continuation::AnasaziGroup::computeEigenvalues:"
	 << endl
         <<  "Anasazi Eigensolver requested but not compiled in!" << endl;
  }
  return NOX::Abstract::Group::Ok;
#endif
}


NOX::Abstract::Group::ReturnType
LOCA::Continuation::AnasaziGroup::applyAnasaziOperator(NOX::Parameter::List &params, 
                           const NOX::Abstract::Vector &input, 
                           NOX::Abstract::Vector &result) const
{

     const LOCA::TimeDependent::AbstractGroup& timeGroup = 
      dynamic_cast<const LOCA::TimeDependent::AbstractGroup&>(*this);

  // Apply Cayley transformation 

  if(cayley)   {

   // First, apply shifted matrix (J - mu M) to input

    const double& shift = -1.00*mu;
 
    NOX::Abstract::Vector *result2 = result.clone(NOX::ShapeCopy);

    timeGroup.applyShiftedMatrix(input,*result2,shift,hasMassMatrix);
    
   // Second, apply inverse of shifted matrix (J - sigma M)
    const double& shift2 = -1.00*sigma;

    timeGroup.applyShiftedMatrixInverse(*result2,result,shift2,hasMassMatrix,params);
    delete result2;
    return NOX::Abstract::Group::Ok;   
}
  else {  

    double omega = -1.00*shift_;

    // Shift and invert if no Cayley transformation

    // If you have a mass matrix, first apply this to input

    if(hasMassMatrix)   {
    NOX::Abstract::Vector *input2 = input.clone(NOX::ShapeCopy);
    
    timeGroup.applyMassMatrix(input, *input2);  
    
    const NOX::Abstract::Vector *massinput = input2->clone(NOX::DeepCopy);
    
   // Second, apply inverse of shifted matrix (J - shift_  M)
    timeGroup.applyShiftedMatrixInverse(*massinput,result,omega,hasMassMatrix,params);

    delete input2;
    delete massinput;

    }
    
    else   {
         timeGroup.applyShiftedMatrixInverse(input,result,omega,hasMassMatrix,params);
    }
  }
}
