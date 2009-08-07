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
#include "GForwardSolve.h"
#include "GAdjointSolve.h"
#include "GErrorEstimate.h"
#include "GMeshRefine.h"
#include "GModel.h"

using namespace GAASP;

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
	double **fSoln;			// Forward solution
	double *mesh;			// Forward solution mesh
	double **fDer;			// Forward solution derivatives
	double **aSoln;			// Adjoint solution
	double **intError;		// Interval contributions
	double errorEstimate;

	//----------------------------------------------------------------------------
	//	Lorenz Test
	//----------------------------------------------------------------------------

	// Define the simulation control values.
	double sTime = 0.0;
	double eTime = 10.0;
	double timeStep = 0.001;
	double uTOL = 0.00001;
	int nSteps;
	TMethodFlag method = Method_DG1;
	TSimControl fsim = {sTime, eTime, timeStep, 0, method};
	TAdjointControl asim = {sTime, eTime, timeStep, 0, Error_Avg, method, 1};

	// Define the model and the initial conditions.
	params.push_back(10.0);
	params.push_back(28.0);
	params.push_back(8.0/3.0);
	shared_ptr<Lorenz_GModel> model(new Lorenz_GModel(params));
	initialConditions.push_back(-6.9742);
	initialConditions.push_back(-7.008);
	initialConditions.push_back(25.1377);

	//	Forward, adjoint, and error estimates on the coarse mesh.
	GForwardSolve forwardSolver(model,fsim,initialConditions);
	forwardSolver.solve();
	mesh = forwardSolver.getMesh();
	forwardSolver.printSolution();
	nSteps = forwardSolver.getNSteps();
	fSoln = forwardSolver.getSolution();
	fDer = forwardSolver.getDerivatives();

	GAdjointSolve adjointSolver(model, asim, fSoln, mesh);
	adjointSolver.solve();
	aSoln = adjointSolver.getSolution();
	adjointSolver.printSolution();

	GErrorEstimate errorEstimator(fSoln, aSoln, fDer, mesh, initialConditions.size(), nSteps, model, method);
	errorEstimator.compute();
	errorEstimator.printErrorEstimates();

	return 0;
}
	

