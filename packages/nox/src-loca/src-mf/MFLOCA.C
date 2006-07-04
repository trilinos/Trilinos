// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
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
#include <MFError.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void MFSetError(int,char*,char*,int,char*);
void MFLOCAFreeData(void *data);

int MFProjectLOCA(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*);
int MFTangentLOCA(int,int,MFNVector,MFNKMatrix,void*);
double MFScaleLOCA(int,int,MFNVector,MFNKMatrix,void*);

int MFLOCAProjectToSave(MFNVector,double*,void*);
int MFLOCAProjectToDraw(MFNVector,double*,void*);
int MFLOCAProjectForBB(MFNVector,double*,void*);

double MFPrintMetricLOCA(double*,double*);

LOCAData::LOCAData(
     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
     const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& top_params,
     const Teuchos::RefCountPtr<NOX::Solver::Generic>& s, 
     const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>& g, 
     const Teuchos::RefCountPtr<Teuchos::ParameterList>& par,
     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& st,
     const Teuchos::RefCountPtr< list<ParamData> >& conParamData) :
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
  Teuchos::RefCountPtr<Teuchos::ParameterList> stepperList = 
    topParams->getSublist("Stepper");
  maxNonlinearIterations = 
    static_cast<double>(stepperList->get("Max Nonlinear Iterations", 
						 15));
  aggressiveness = stepperList->get("Aggressiveness", 0.0);
  solutionMax = stepperList->get("Max Solution Component", 1.0e16);
}

MFImplicitMF MFIMFCreateLOCA(LOCAData* data)
 {
  MFImplicitMF loca;
  MFNSpace space;

  loca=MFIMFCreateBaseClass(-1, data->np, "LOCA");

  space=MFCreateLOCANSpace(data);
  MFIMFSetSpace(loca,space);
  MFFreeNSpace(space);

  MFIMFSetData(loca,(void*)data);
  data->space=space;
  MFRefNSpace(space);
  MFIMFSetFreeData(loca,MFLOCAFreeData);
  MFIMFSetProject(loca,MFProjectLOCA);
  MFIMFSetTangent(loca,MFTangentLOCA);
  MFIMFSetScale(loca,MFScaleLOCA);
  MFIMFSetProjectForSave(loca,MFLOCAProjectToDraw);
  MFIMFSetProjectForDraw(loca,MFLOCAProjectToDraw);
  MFIMFSetProjectForBB(loca,MFLOCAProjectToDraw);

  return loca;
 }

void MFLOCAFreeData(void *d)
 {
    LOCAData* data = (LOCAData *)d;
    delete data;
 }

int MFLOCAProjectToDraw(MFNVector u, double *x, void *d)
 {
  LOCAData* data = (LOCAData *)d; 

  if(x==(double*)NULL)
    return data->grp->projectToDrawDimension();

  LOCANVectorData* v_data = (LOCANVectorData *)MFNVectorGetData(u);

  data->grp->projectToDraw(*(v_data->u_ptr), x);

  return 0;
 }

int MFProjectLOCA(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,
		  void *d,int *index)
{
  static int stepNumber = 1;
  int i;
  *index = 1;

  LOCAData* data = (LOCAData *)d;
  for (i=0; i<k; i++) {
    MFNVector tmp =  MFMColumn(mPhi,i);
    LOCANVectorData* tmp2_data = (LOCANVectorData *) MFNVectorGetData(tmp);
    data->grp->setPredictorTangentDirection(*(tmp2_data->u_ptr), i);
    MFFreeNVector(tmp);
  }
  
  LOCANVectorData* u0_data = (LOCANVectorData *) MFNVectorGetData(vu0);
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
    list<ParamData>::iterator it = data->paramData->begin();
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
  data->solver->reset(data->grp, data->status, 
		      Teuchos::rcp(&(data->p->sublist("NOX")),false));
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
      list<ParamData>::iterator it = data->paramData->begin();
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
    LOCANVectorData* u_data = (LOCANVectorData *) MFNVectorGetData(vu);
     
    *(Teuchos::rcp_dynamic_cast<NOX::Abstract::Group>(data->grp)) = 
      data->solver->getSolutionGroup();
    *(u_data->u_ptr) = data->grp->getX(); /* overloaded deep copy */
    data->grp->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
    
    if (data->globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      data->globalData->locaUtils->out() 
	<< "\n" << data->globalData->locaUtils->fill(72, '~') << "\n";
      data->globalData->locaUtils->out() 
	<< "End of Continuation Step " << stepNumber << " : " << std::endl;
      list<ParamData>::iterator it = data->paramData->begin();
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


int MFTangentLOCA(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d)
{
   LOCAData* data = (LOCAData *)d;

   LOCANVectorData* u0_data = (LOCANVectorData *) MFNVectorGetData(vu);
   data->grp->setX(*(u0_data->u_ptr));
   data->grp->computePredictor();

   const LOCA::MultiContinuation::ExtendedMultiVector& pred = 
     data->grp->getPredictorTangent();

   for (int i=0; i<k; i++) {
     Teuchos::RefCountPtr<LMCEV> t = 
       Teuchos::rcp_dynamic_cast<LMCEV>(pred[i].clone());
     MFNVector tmp =  MFCreateLOCANVectorWithData(t);
     MFMSetColumn(mPhi, i, tmp);
     MFFreeNVector(tmp);
   }

   MFGramSchmidt(data->space,mPhi);

   return 1;
}

double MFScaleLOCA(int n,int k,MFNVector u,MFNKMatrix Phi,void *d)
{
  LOCAData* data = (LOCAData *)d;
  if (data->radius < 0.0) {
    data->radius = 0.0;
    data->minRadius = 0.0;
    data->maxRadius = 0.0;
    list<ParamData>::iterator it = data->paramData->begin();
    for (int i=0; i<k; i++) {
      double dpidsj = 0.0;
      for (int j=0; j<k; j++) {
	MFNVector tmp =  MFMColumn(Phi,j);
	LOCANVectorData* tmp2_data = (LOCANVectorData *) MFNVectorGetData(tmp);
	dpidsj += fabs(tmp2_data->u_ptr->getScalar(i));
	MFFreeNVector(tmp);
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
