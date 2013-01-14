// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "MFLOCA.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"

extern "C" {

#include <MFPrint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void MFLOCAFreeData(void *data,MFErrorHandler);

int MFProjectLOCA(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,
		  MFErrorHandler);
int MFTangentLOCA(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
double MFScaleLOCA(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);

int MFLOCAProjectToSave(MFNVector,double*,void*,MFErrorHandler);
int MFLOCAProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
int MFLOCAProjectForBB(MFNVector,double*,void*,MFErrorHandler);

double MFPrintMetricLOCA(double*,double*,MFErrorHandler);

LOCAData::LOCAData(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& top_params,
     const Teuchos::RCP<NOX::Solver::Generic>& s, 
     const Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>& g, 
     const Teuchos::RCP<Teuchos::ParameterList>& par,
     const Teuchos::RCP<NOX::StatusTest::Generic>& st,
     const Teuchos::RCP< std::list<ParamData> >& conParamData) :
  globalData(global_data),
  topParams(top_params),
  solver(s), 
  grp(g), 
  p(par), 
  status(st),
  paramData(conParamData),
  space(NULL), 
  np(g->getNumParams()),
  maxNonlinearIterations(1.0),
  aggressiveness(0.0),
  radius(-1.0),
  maxRadius(0.0),
  minRadius(0.0),
  solutionMax(0.0)
{
  Teuchos::RCP<Teuchos::ParameterList> stepperList = 
    topParams->getSublist("Stepper");
  maxNonlinearIterations = 
    static_cast<double>(stepperList->get("Max Nonlinear Iterations", 
						 15));
  aggressiveness = stepperList->get("Aggressiveness", 0.0);
  solutionMax = stepperList->get("Max Solution Component", 1.0e16);
  mfErrorHandler = MFCreateErrorHandler();
}

MFImplicitMF MFIMFCreateLOCA(LOCAData* data)
 {
  MFImplicitMF loca;
  MFNSpace space;

  loca=MFIMFCreateBaseClass(-1, data->np, "LOCA", data->mfErrorHandler);

  space=MFCreateLOCANSpace(data);
  MFIMFSetSpace(loca,space, data->mfErrorHandler);
  MFFreeNSpace(space, data->mfErrorHandler);

  MFIMFSetData(loca,(void*)data, data->mfErrorHandler);
  data->space=space;
  MFRefNSpace(space, data->mfErrorHandler);
  MFIMFSetFreeData(loca,MFLOCAFreeData, data->mfErrorHandler);
  MFIMFSetProject(loca,MFProjectLOCA, data->mfErrorHandler);
  MFIMFSetTangent(loca,MFTangentLOCA, data->mfErrorHandler);
  MFIMFSetScale(loca,MFScaleLOCA, data->mfErrorHandler);
  MFIMFSetProjectForSave(loca,MFLOCAProjectToDraw, data->mfErrorHandler);
  MFIMFSetProjectForDraw(loca,MFLOCAProjectToDraw, data->mfErrorHandler);
  MFIMFSetProjectForBB(loca,MFLOCAProjectToDraw, data->mfErrorHandler);

  return loca;
 }

void MFLOCAFreeData(void *d, MFErrorHandler err)
 {
    LOCAData* data = (LOCAData *)d;
    delete data;
 }

int MFLOCAProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler err)
 {
  LOCAData* data = (LOCAData *)d; 

  if(x==(double*)NULL)
    return data->grp->projectToDrawDimension();

  LOCANVectorData* v_data = (LOCANVectorData *)MFNVectorGetData(u, err);

  data->grp->projectToDraw(*(v_data->u_ptr), x);

  return 0;
 }

int MFProjectLOCA(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,
		  void *d,int *index, MFErrorHandler err)
{
  static int stepNumber = 1;
  int i;
  *index = 1;

  LOCAData* data = (LOCAData *)d;
  for (i=0; i<k; i++) {
    MFNVector tmp =  MFMColumn(mPhi,i, err);
    LOCANVectorData* tmp2_data = (LOCANVectorData *) MFNVectorGetData(tmp, err);
    data->grp->setPredictorTangentDirection(*(tmp2_data->u_ptr), i);
    MFFreeNVector(tmp, err);
  }
  
  LOCANVectorData* u0_data = (LOCANVectorData *) MFNVectorGetData(vu0, err);
  data->grp->setPrevX(*(u0_data->u_ptr));
  data->grp->setX(*(u0_data->u_ptr));
  for (i=0; i<k; i++) {
    data->grp->setStepSize(0.0, i);
  }

  if (data->globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    data->globalData->locaUtils->out() 
      << "\n" << data->globalData->locaUtils->fill(72, '~') << "\n";
    data->globalData->locaUtils->out() 
      << "Start of Continuation Step " << stepNumber <<" : " << std::endl;
    std::list<ParamData>::iterator it = data->paramData->begin();
    for (i=0; i<k; i++) {
      data->globalData->locaUtils->out() 
	<< "\tParameter: " << it->name << " = " 
	<< data->globalData->locaUtils->sciformat(data->grp->getContinuationParameter(i))
	<< std::endl;
      it++;
    }
    data->globalData->locaUtils->out() 
      << data->globalData->locaUtils->fill(72, '~') << "\n" << std::endl;
  }

  data->grp->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

  data->grp->computeF();
  data->solver->reset(data->grp->getX());
  NOX::StatusTest::StatusType status = data->solver->solve();
    
  if (status != NOX::StatusTest::Converged) {
    if (data->globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      data->globalData->locaUtils->out() 
	<< std::endl << data->globalData->locaUtils->fill(72, '~') 
	<< std::endl;
      data->globalData->locaUtils->out() 
	<< "Continuation Step Number " << stepNumber 
	<< " experienced a convergence failure in\n"
	<< "the nonlinear solver after "<< data->solver->getNumIterations() 
	<<" Iterations\n";
      data->globalData->locaUtils->out() 
	<< "Value of continuation parameters at failed step:" << std::endl;
      std::list<ParamData>::iterator it = data->paramData->begin();
      for (i=0; i<k; i++) {
	data->globalData->locaUtils->out() 
	  << "\tParameter: " << it->name << " = " 
	  << data->globalData->locaUtils->sciformat(data->grp->getContinuationParameter(i))
	  << std::endl;
	it++;
      }
      data->globalData->locaUtils->out() 
	<< data->globalData->locaUtils->fill(72, '~') << std::endl;
    }
    data->grp->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);
    return 0;
  }
  else {
    LOCANVectorData* u_data = (LOCANVectorData *) MFNVectorGetData(vu, err);
     
    *(Teuchos::rcp_dynamic_cast<NOX::Abstract::Group>(data->grp)) = 
      data->solver->getSolutionGroup();
    *(u_data->u_ptr) = data->grp->getX(); /* overloaded deep copy */
    data->grp->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
    
    if (data->globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      data->globalData->locaUtils->out() 
	<< "\n" << data->globalData->locaUtils->fill(72, '~') << "\n";
      data->globalData->locaUtils->out() 
	<< "End of Continuation Step " << stepNumber << " : " << std::endl;
      std::list<ParamData>::iterator it = data->paramData->begin();
      for (i=0; i<k; i++) {
	data->globalData->locaUtils->out() 
	  << "\tParameter: " << it->name << " = " 
	  << data->globalData->locaUtils->sciformat(data->grp->getContinuationParameter(i))
	  << std::endl;
	it++;
      }
      data->globalData->locaUtils->out() 
	<< "--> Step Converged in "
	<< data->solver->getNumIterations() 
	<<" Nonlinear Solver Iterations!\n";
      data->globalData->locaUtils->out() 
	<< data->globalData->locaUtils->fill(72, '~') << "\n" << std::endl;
    }
    ++stepNumber;
    return 1;
  }
}


int MFTangentLOCA(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, 
		  MFErrorHandler err)
{
   LOCAData* data = (LOCAData *)d;

   LOCANVectorData* u0_data = (LOCANVectorData *) MFNVectorGetData(vu, err);
   data->grp->setX(*(u0_data->u_ptr));
   data->grp->computePredictor();

   const LOCA::MultiContinuation::ExtendedMultiVector& pred = 
     data->grp->getPredictorTangent();

   for (int i=0; i<k; i++) {
     Teuchos::RCP<LMCEV> t = 
       Teuchos::rcp_dynamic_cast<LMCEV>(pred[i].clone());
     MFNVector tmp =  MFCreateLOCANVectorWithData(t,err);
     MFMSetColumn(mPhi, i, tmp, err);
     MFFreeNVector(tmp, err);
   }

   MFGramSchmidt(data->space,mPhi, err);

   return 1;
}

double MFScaleLOCA(int n,int k,MFNVector u,MFNKMatrix Phi,void *d, 
		   MFErrorHandler err)
{
  LOCAData* data = (LOCAData *)d;
  if (data->radius < 0.0) {
    data->radius = 0.0;
    data->minRadius = 0.0;
    data->maxRadius = 0.0;
    std::list<ParamData>::iterator it = data->paramData->begin();
    for (int i=0; i<k; i++) {
      double dpidsj = 0.0;
      for (int j=0; j<k; j++) {
	MFNVector tmp =  MFMColumn(Phi,j, err);
	LOCANVectorData* tmp2_data = (LOCANVectorData *) MFNVectorGetData(tmp, err);
	dpidsj += fabs(tmp2_data->u_ptr->getScalar(i));
	MFFreeNVector(tmp, err);
      }
      data->radius += it->initialStepSize / dpidsj;
      data->minRadius += it->minStepSize / dpidsj;
      data->maxRadius += it->maxStepSize / dpidsj;
      it++;
    }
    data->radius /= (double) k;
    data->minRadius /= (double) k;
    data->maxRadius /= (double) k;
  }
  else {
    NOX::StatusTest::StatusType status = data->solver->getStatus();
    if (status != NOX::StatusTest::Converged) {
      data->radius *= 0.7;
      if (data->radius < data->minRadius)
	data->globalData->locaErrorCheck->throwError(
						   "MFScaleLOCA()",
						   "Reached minimum radius!");
    }
    else {
      double numNonlinearSteps = 
	static_cast<double>(data->solver->getNumIterations());
      double factor = (data->maxNonlinearIterations - numNonlinearSteps) 
	/ (data->maxNonlinearIterations);
      data->radius *= (1.0 + data->aggressiveness * factor * factor);
      if (data->radius > data->maxRadius)
	data->radius = data->maxRadius;
    }
  }   

  data->globalData->locaUtils->out() 
    << "radius = " << data->radius << std::endl;
  return data->radius; 
}

} /* extern C */
