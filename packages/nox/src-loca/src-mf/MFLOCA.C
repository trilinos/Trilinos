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

#include <MFLOCA.H>
#include "LOCA_Utils.H"

extern "C" {

#include <MFPrint.h>
#include <MFError.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static char MFLOCAMFErrorMsg[256]="";
void MFSetError(int,char*,char*,int,char*);
void MFLOCAFreeData(void *data);

int MFProjectLOCA(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*);
int MFTangentLOCA(int,int,MFNVector,MFNKMatrix,void*);
double MFScaleLOCA(int,int,MFNVector,MFNKMatrix,void*);

int MFLOCAProjectToSave(MFNVector,double*,void*);
int MFLOCAProjectToDraw(MFNVector,double*,void*);
int MFLOCAProjectForBB(MFNVector,double*,void*);

double MFPrintMetricLOCA(double*,double*);

LOCAData::LOCAData(NOX::Solver::Generic& s, 
		   LOCA::MultiContinuation::ExtendedGroup& g, 
		   NOX::Parameter::List& par,
		   NOX::StatusTest::Generic& st,
		   list<ParamData>& conParamData) :
  solver(s), 
  grp(g), 
  p(par), 
  status(st),
  paramData(conParamData),
  space(NULL), 
  np(g.getNumParams()),
  maxNonlinearIterations(1.0),
  aggressiveness(0.0),
  radius(-1.0),
  maxRadius(0.0),
  minRadius(0.0),
  solutionMax(0.0)
{
  NOX::Parameter::List& stepperList = LOCA::Utils::getSublist("Stepper");
  maxNonlinearIterations = 
    static_cast<double>(stepperList.getParameter("Max Nonlinear Iterations", 
						 15));
  aggressiveness = stepperList.getParameter("Aggressiveness", 0.0);
  solutionMax = stepperList.getParameter("Max Solution Component", 1.0e16);
}

MFImplicitMF MFIMFCreateLOCA(LOCAData* data)
 {
  static char RoutineName[]={"MFIMFCreateLOCA"};
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
    return data->grp.projectToDrawDimension();

  LMCEV* v = (LMCEV *)MFNVectorGetData(u);

  data->grp.projectToDraw(*v, x);

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
    LMCEV* tmp2 = (LMCEV *) MFNVectorGetData(tmp);
    data->grp.setPredictorDirection(*tmp2, i);
    MFFreeNVector(tmp);
  }
  
  LMCEV* u0 = (LMCEV *) MFNVectorGetData(vu0);
  data->grp.setPrevX(*u0);
  data->grp.setX(*u0);
  for (i=0; i<k; i++) {
    data->grp.setStepSize(0.0, i);
  }

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
    cout << "\n" << LOCA::Utils::fill(72, '~') << "\n";
    cout << "Start of Continuation Step " << stepNumber <<" : " << endl;
    list<ParamData>::iterator it = data->paramData.begin();
    for (i=0; i<k; i++) {
      cout << "\tParameter: " << it->name << " = " 
	   << LOCA::Utils::sci(data->grp.getContinuationParameter(i))
	   << endl;
      it++;
    }
    cout << LOCA::Utils::fill(72, '~') << "\n" << endl;
  }

  data->grp.computeF();
  data->solver.reset(data->grp, data->status, data->p.sublist("NOX"));
  NOX::StatusTest::StatusType status = data->solver.solve();
    
  if (status != NOX::StatusTest::Converged) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
      cout << endl << LOCA::Utils::fill(72, '~') << endl;
      cout << "Continuation Step Number " << stepNumber 
           << " experienced a convergence failure in\n"
           << "the nonlinear solver after "<< data->solver.getNumIterations() 
	   <<" Iterations\n";
      cout << "Value of continuation parameters at failed step:" << endl;
      list<ParamData>::iterator it = data->paramData.begin();
      for (i=0; i<k; i++) {
	cout << "\tParameter: " << it->name << " = " 
	     << LOCA::Utils::sci(data->grp.getContinuationParameter(i))
	     << endl;
	it++;
      }
      cout << LOCA::Utils::fill(72, '~') << endl;
    }
    return 0;
  }
  else {
    LMCEV* u = (LMCEV *) MFNVectorGetData(vu);
     
    data->grp = data->solver.getSolutionGroup();
    *u = data->grp.getX(); /* overloaded deep copy */
    data->grp.notifyCompletedStep();
    
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
      cout << "\n" << LOCA::Utils::fill(72, '~') << "\n";
      cout << "End of Continuation Step " << stepNumber << " : " << endl;
      list<ParamData>::iterator it = data->paramData.begin();
      for (i=0; i<k; i++) {
	cout << "\tParameter: " << it->name << " = " 
	     << LOCA::Utils::sci(data->grp.getContinuationParameter(i))
	     << endl;
	it++;
      }
      cout << "--> Step Converged in "
           << data->solver.getNumIterations() 
	   <<" Nonlinear Solver Iterations!\n";
      cout << LOCA::Utils::fill(72, '~') << "\n" << endl;

      //u->print();
    }
    ++stepNumber;
    return 1;
  }
}


int MFTangentLOCA(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d)
{
   LOCAData* data = (LOCAData *)d;

   LMCEV* u0 = (LMCEV *) MFNVectorGetData(vu);
   data->grp.setX(*u0);
   data->grp.computePredictor();

   const LOCA::MultiContinuation::ExtendedMultiVector& pred = 
     data->grp.getPredictorDirections();

   for (int i=0; i<k; i++) {
     LMCEV* t = dynamic_cast<LMCEV*>(pred[i].clone());
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
    list<ParamData>::iterator it = data->paramData.begin();
    for (int i=0; i<k; i++) {
      double dpidsj = 0.0;
      for (int j=0; j<k; j++) {
	MFNVector tmp =  MFMColumn(Phi,j);
	LMCEV* tmp2 = (LMCEV *) MFNVectorGetData(tmp);
	dpidsj += fabs(tmp2->getScalar(i));
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
    NOX::StatusTest::StatusType status = data->solver.getStatus();
    if (status != NOX::StatusTest::Converged) {
      data->radius *= 0.7;
      if (data->radius < data->minRadius)
	LOCA::ErrorCheck::throwError("MFScaleLOCA()",
				     "Reached minimum radius!");
    }
    else {
      double numNonlinearSteps = 
	static_cast<double>(data->solver.getNumIterations());
      double factor = (data->maxNonlinearIterations - numNonlinearSteps) 
	/ (data->maxNonlinearIterations);
      data->radius *= (1.0 + data->aggressiveness * factor * factor);
      if (data->radius > data->maxRadius)
	data->radius = data->maxRadius;
    }
  }   

  cout << "radius = " << data->radius << endl;
  return data->radius; 
}

} /* extern C */
