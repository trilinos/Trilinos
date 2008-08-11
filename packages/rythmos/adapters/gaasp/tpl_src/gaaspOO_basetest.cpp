// ----------------------------------------------------------------------------
//	Developed with funding from the National Science Foundation
//	grant NSF EF-0434354
// ----------------------------------------------------------------------------
//	Organization:	Department of Mathematics
//			Colorado State University, Fort Collins, CO 80523 USA
//	Project:  N.S.F. Mathematical and Statistical Tools for Understanding
//	          Complex Systems in the Environment
//	File:	  gaaspOO.cpp
//	Class:	  
//
//	Description:
//	Example main program for using the GForwardSolve class.  
// ----------------------------------------------------------------------------
//	Author:	Jeff Sandelin, sandelin@math.colostate.edu, 17Apr2008
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cstdlib>
#include <ctime>
#include <climits>

//#include "initGaasp.h"
#include "MemoryMgmt.h"
#include "gforwardsolve.h"
#include "gadjointsolve.h"
#include "gerrorestimate.h"
#include "gmeshrefine.h"
#include "GModel.h"

namespace GAASP {

using std::cout;
using std::endl;

// Function parameter list and initial conditions
//  NOTE: A slight difference in how these are handled. The functions
//        can access these directly, but the initial conditions are passed
//        to the solver since you have to call the function with the current
//        state for the timestep you are at. In some cases, the initial  be
//        conditions may be used in the true solution function as a fixed
//        parameter.
std::vector<double> params;
std::vector<double> initialConditions;

// Utility functions or prototypes
void DisplayArgs (
	int const argc,
	char const * const * const argv)
{
	if ( argc > 1 )
	{
	    cout << "Args: ";
	    for ( int i = 1; i < argc; ++i )
		cout << argv[i] << " ";
	    cout << endl;
	}
}

//---	Main
int main(int argc, char *argv[])
{
	
	// Define the simulation control values.
	double startTime = 0.0;
	double endTime = 10.0;
	double timeStep = 0.01;
	double uTOL = 0.00001;
	int nsteps = 0;
	TMethodFlag method = Method_DG1;
	TSimControl fsim = {startTime, endTime, timeStep, nsteps, method};
	TAdjointControl asim = {startTime, endTime, timeStep, nsteps, Error_Avg, Method_DG1, 1};

	// Define the model and the initial conditions.
	GModel model;
	initialConditions.push_back(-6.9742);
	initialConditions.push_back(-7.008);
	initialConditions.push_back(25.1377);
	int dim = initialConditions.size();

	// Define the return variables.
	double **fsoln;		// Forward solution
	double *mesh;		// Forward solution mesh
	double **sder;		// Forward solution derivatives
	double **asoln;		// Adjoint solution
	double errEstimate;	// Error estimate
	double **intError;	// Interval contributions
	double *newMesh;	// Refined mesh
	int nNodes;			// Number of nodes in refined mesh

	// Initial class instantiations and solve.
	GForwardSolve forwardSolver(model,fsim,initialConditions,params);
	forwardSolver.solve();
	mesh = forwardSolver.getMesh();
	nsteps = forwardSolver.getNSteps();
	fsoln = forwardSolver.getSolution();
	sder = forwardSolver.getDerivatives();

	GAdjointSolve adjointSolver(model, asim, fsoln, mesh);
	adjointSolver.solve();
	asoln = adjointSolver.getSolution();

	GErrorEstimate errorEstimator(fsoln, asoln, sder, mesh, initialConditions.size(), nsteps, model, method);
	errorEstimator.compute();
	errEstimate = errorEstimator.getErrorEstimate();
	intError = errorEstimator.getIntervalContributions();

	GMeshRefine meshRefiner(mesh, intError, Refine_WeightedAdaptive, Method_DG1, uTOL, dim, nsteps+1);

	while(fabs(errEstimate) > uTOL)
	{
		meshRefiner.refine();
		newMesh = meshRefiner.getNewMesh();
		nNodes = meshRefiner.getNumberofNodes();
		forwardSolver.updateMesh(newMesh,nNodes);
		forwardSolver.solve();
		forwardSolver.printSolution();
		mesh = forwardSolver.getMesh();
		nsteps = forwardSolver.getNSteps();
		fsoln = forwardSolver.getSolution();
		sder = forwardSolver.getDerivatives();
		adjointSolver.updateData(mesh,nsteps,fsoln);
		adjointSolver.solve();
		adjointSolver.printSolution();
		asoln = adjointSolver.getSolution();
		errorEstimator.updateData(mesh,nsteps,fsoln,sder,asoln);
		errorEstimator.compute();
		errEstimate = errorEstimator.getErrorEstimate();
		break;
	}




	/*
	// Hard coded for a maximum of 1 iterations.
	int refinementCycle = 0;
	int maxRefinementCycles = 1;
	for(;;)
	{

		// Solve the forward problem
		forwardSolver.solve();
		fsoln = forwardSolver.getSolution();
		sder = forwardSolver.getDerivatives();
		nsteps = forwardSolver.getNSteps();

		// Solve the adjoint problem
		adjointSolver.updateData(mesh,nsteps+1,fsoln);
		adjointSolver.solve();
		asoln = adjointSolver.getSolution();

		// Compute the error estimate.
		errorEstimator.updateData(mesh, nsteps+1, fsoln, sder, asoln);
		errorEstimator.compute();
		errEstimate = errorEstimator.getErrorEstimate();
		if(fabs(errEstimate) < uTOL || refinementCycle == maxRefinementCycles)
		{
			forwardSolver.printSolution();
			break;
		}
		refinementCycle++;
		intError = errorEstimator.getIntervalContributions();

		// Mesh refinement
		
		meshRefiner.refine();

		// Do the forward solve using the new mesh.
		newMesh = meshRefiner.getNewMesh();
		nNodes = meshRefiner.getNumberofNodes();
		forwardSolver.updateMesh(newMesh,nNodes);

	}
	*/
    return 0;
}
	
} // namespace GAASP

