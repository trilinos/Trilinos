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
	// Define the return variables for coarse mesh.
	double **coarseFSoln;			// Forward solution
	double *coarseMesh;				// Forward solution mesh
	double **coarseDer;				// Forward solution derivatives
	double **coarseASoln;			// Adjoint solution
	double **coarseIntError;		// Interval contributions
	double errorEstimate;

	// Define the return variables for partitions.
	double **partFSoln;				// Forward solution
	double *partMesh;				// Forward solution mesh
	double **partDer;				// Forward solution derivatives
	double **partASoln;				// Adjoint solution
	double **partIntError;			// Interval contributions
	int nPartSteps;

	// Define the return variables for fine mesh.
	double **fineFSoln;				// Forward solution
	double *fineMesh;				// Forward solution mesh
	double **fineDer;				// Forward solution derivatives
	double **fineASoln;				// Adjoint solution
	double **fineIntError;			// Interval contributions
	int nFineSteps;

	//----------------------------------------------------------------------------
	//	Test for domain decomposition.
	//----------------------------------------------------------------------------

	// Define the simulation control values.
	double sTime = 0.0;
	double eTime = 10.0;
	double coarseTimeStep = 0.1;
	double uTOL = 0.00001;
	int nCoarseSteps = 0;
	TMethodFlag method = Method_DG1;
	TSimControl fsim = {sTime, eTime, coarseTimeStep, 0, method};
	TAdjointControl asim = {sTime, eTime, coarseTimeStep, 0, Error_Avg, Method_DG1, 1};

	// Define the model and the initial conditions.
	params.push_back(10.0);
	params.push_back(28.0);
	params.push_back(8.0/3.0);
	GModel model(params);
	initialConditions.push_back(-6.9742);
	initialConditions.push_back(-7.008);
	initialConditions.push_back(25.1377);
	int dim = initialConditions.size();

	//	Forward, adjoint, and error estimates on the coarse mesh.
	GForwardSolve forwardSolver(model,fsim,initialConditions);
	forwardSolver.solve();
	coarseMesh = forwardSolver.getMesh();
	forwardSolver.printSolution();
	nCoarseSteps = forwardSolver.getNSteps();
	coarseFSoln = forwardSolver.getSolution();
	coarseDer = forwardSolver.getDerivatives();

	GAdjointSolve adjointSolver(model, asim, coarseFSoln, coarseMesh);
	adjointSolver.solve();
	coarseASoln = adjointSolver.getSolution();

	GErrorEstimate errorEstimator(coarseFSoln, coarseASoln, coarseDer, coarseMesh, initialConditions.size(), nCoarseSteps, model, method);
	errorEstimator.compute();
	coarseIntError = errorEstimator.getIntervalContributions();

	//	Fine mesh info that doesn't change.
	double timeStep = 0.01;
	double factor = 10;
	int numPartitions = 2;

	//	Compute the start and end times for the partitions.
	double *startTime = 0;
	double *endTime = 0;
	startTime = CreateArray<double>(numPartitions);
	endTime = CreateArray<double>(numPartitions);
	for(int i=0; i<numPartitions; i++)
	{
		startTime[i] = coarseMesh[int(nCoarseSteps/numPartitions)*i];
		endTime[i] = coarseMesh[int(nCoarseSteps/numPartitions)*(i+1)];
	}
	nFineSteps = (endTime[numPartitions-1]-startTime[0])/timeStep;

	//	Compute initial conditions for the partitions.
	double **IC;
	IC = CreateArray<double>(numPartitions,dim);
	for(int i=0; i<numPartitions; i++)
	{
		for(int j=0; j<dim; j++)
		{
			IC[i][j] = coarseFSoln[j][int(nCoarseSteps/numPartitions)*i];
		}
	}

	//	Loop through the partitions.
	for(int i=0; i<numPartitions; i++)
	{
		//	Set simulation control for partition.
		nPartSteps = (endTime[i]-startTime[i])/timeStep;

		// Initial conditions.
		initialConditions.clear();
		initialConditions.push_back(IC[i][0]);
		initialConditions.push_back(IC[i][1]);
		initialConditions.push_back(IC[i][2]);
		dim = initialConditions.size();

		//	Generate the partition mesh.
		partMesh = CreateArray<double>(nCoarseSteps*factor/numPartitions+1);
		for(int j=0; j<nPartSteps+1; j++)
		{
			partMesh[j] = startTime[i] + (j*timeStep);
		}

		//	Update data and do the solves.
		forwardSolver.updateIC(initialConditions);
		forwardSolver.updateMesh(partMesh,nPartSteps+1);
		forwardSolver.solve();
		forwardSolver.printSolution(i+1);
		partFSoln = forwardSolver.getSolution();
		partDer = forwardSolver.getDerivatives();
		adjointSolver.updateData(partMesh,nPartSteps,partFSoln);
		adjointSolver.solve();
		adjointSolver.printSolution(i+1);
		partASoln = adjointSolver.getSolution();
		errorEstimator.updateData(partMesh,nPartSteps,partFSoln,partDer,partASoln);
		errorEstimator.compute();
		errorEstimate = errorEstimator.getErrorEstimate();
	}
	return 0;
}
	
} // namespace GAASP

