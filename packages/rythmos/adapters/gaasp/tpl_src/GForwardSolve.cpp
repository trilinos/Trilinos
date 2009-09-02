// ----------------------------------------------------------------------------
//	Developed with funding from the Department of Energy
//  Grants DE-FG02-04ER25620, DE-FG02-05ER25699, DE-FC02-07ER54909
//	and
//	National Science Foundation
//	Grants DMS-0715135, DGE-0221595003, MSPA-CSE-0434354
// ----------------------------------------------------------------------------
//	Copyright (c) 2007 Colorado State University. All rights reserved.
// ----------------------------------------------------------------------------
//	Organization:	Department of Mathematics - Estep Research Group
//					Colorado State University, Fort Collins, CO 80523 USA
//					www.math.colostate.edu/~estep/gaasp
//	Project:  Globally Accurate Adaptive Sensitivity Package
//	File:	  GForwardSolve.cpp
//	Class:	  GForwardSolve
//
//	Description:
//	Handles solving the forward problem for GAASP.
// ----------------------------------------------------------------------------
//	Author:	Jeff Sandelin, sandelin@math.colostate.edu, November, 2007
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------

#include <limits>
#include "assert.h"
#include <string>
#include "GForwardSolve.h"
#include <stdexcept>

namespace GAASP {

using std::numeric_limits;

// ----------------------------------------------------------------------------
//	member constants
// ----------------------------------------------------------------------------

char const * const GForwardSolve::version =		// class version number
	"1.0.0.0";	// major.minor.release.build

typedef void (*TFunction) ( double *, double *, double );

// ----------------------------------------------------------------------------
//	 constructors and destructor
// ----------------------------------------------------------------------------

GForwardSolve::GForwardSolve ()
{
  isSetup_ = false;
}

//---- Model class, TSimControl, initial conditions
GForwardSolve::GForwardSolve (shared_ptr<GModelBase> model, TSimControl simCon, std::vector<double> InitCond)
{
  Setup(model, simCon, InitCond);
}

void GForwardSolve::Setup (shared_ptr<GModelBase> model, TSimControl simCon, std::vector<double> InitCond)
{
	int i;
	modelPtr = model;
	nsteps = int((simCon.endTime-simCon.startTime)/simCon.stepSize);
	sysDim = modelPtr->getDim();
  isSetup_ = true;
	Initialize ();
	tNodes = CreateArray<double>(nsteps+1);
	for( i = 0; i < nsteps+1; i++ )
	{
		tNodes[i] = simCon.startTime + i*simCon.stepSize;
	}
    icCnt = InitCond.size();
	IC = CreateArray<double>(icCnt);
    for ( i = 0; i < icCnt; ++i )
    {
		IC[i] = InitCond[i];
        solution[i][0] = InitCond[i];		// Initial state
    }
	MethodFlag = simCon.methodFlag;
}

//---- Model class, TSimControl, mesh, initial conditions, parameters
GForwardSolve::GForwardSolve (boost::shared_ptr<GModelBase> model, TSimControl simCon, double *mesh, int steps, std::vector<double> IC, std::vector<double> params)
{
  Setup(model, simCon, mesh, steps, IC, params);
}
void GForwardSolve::Setup (boost::shared_ptr<GModelBase> model, TSimControl simCon, double *mesh, int steps, std::vector<double> IC, std::vector<double> params)
{
	int i;
	modelPtr = model;
	nsteps = steps;
	sysDim = modelPtr->getDim();
  isSetup_ = true;
	Initialize ();
	tNodes = CreateArray<double>(nsteps+1);
	for( i = 0; i < nsteps+1; i++ )
	{
		tNodes[i] = mesh[i];
	}
	simCon.startTime = tNodes[0];
	simCon.endTime = tNodes[nsteps];
    int icCnt = IC.size();
    for ( i = 0; i < icCnt; ++i )
    {
        solution[i][0] = IC[i];		// Initial state
    }
	MethodFlag = simCon.methodFlag;
}

GForwardSolve::GForwardSolve (boost::shared_ptr<GModelBase> model, TSimControl simCon, std::vector<double> InitCond, std::vector<double> mesh)
{
  Setup(model,simCon,InitCond,mesh);
}
void GForwardSolve::Setup (boost::shared_ptr<GModelBase> model, TSimControl simCon, std::vector<double> InitCond, std::vector<double> mesh)
{
	modelPtr = model;
	nsteps = mesh.size() - 1;
	sysDim = modelPtr->getDim();
  isSetup_ = true;
	Initialize ();
	tNodes = CreateArray<double>(mesh.size());
	for( unsigned int i = 0; i < mesh.size(); i++ )
	{
		tNodes[i] = mesh[i];
	}
    icCnt = InitCond.size();
	IC = CreateArray<double>(icCnt);
    for ( unsigned int i = 0; i < icCnt; ++i )
    {
		IC[i] = InitCond[i];
        solution[i][0] = InitCond[i];		// Initial state
    }
	MethodFlag = simCon.methodFlag;
}

GForwardSolve::~GForwardSolve ()
{
	Clear ();
}

// ----------------------------------------------------------------------------
//	public functions
// ----------------------------------------------------------------------------

//	Clear
// 	"clear" data members
void GForwardSolve::Clear ()
{
  assert(isSetup_);
	DeleteArray(solution, sysDim);
}

// ----------------------------------------------------------------------------
//	protected functions
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//	private functions
// ----------------------------------------------------------------------------

//	Initialize
void GForwardSolve::Initialize ()
{
  assert(isSetup_);
	machEpsilon = numeric_limits<double>::epsilon();
	allocateSolutionSpace();
}

//	Copy
// 	copy to this
void GForwardSolve::Copy (
	GForwardSolve const & object)
{
	if ( &object )
	{
	}
}

void GForwardSolve::allocateSolutionSpace()
{
  assert(isSetup_);
	solution = CreateArray<double>(sysDim,nsteps+1);
	derivatives = CreateArray<double>(sysDim,nsteps+1);
}

int GForwardSolve::solve()
{
  assert(isSetup_);
//	double const * const t       replace with tNodes
//	int const nsteps             set with size.tNodes
//	int const sysDim					member variable
//  std::vector<double> initialState	member variable
//  initial conditions set in initialize member function

	int n, i;	// Loop counters

	//	Initialize arrays.

    double *x = CreateArray<double>(sysDim);
    double *xout = CreateArray<double>(sysDim);
    double *xp = CreateArray<double>(sysDim * ( MethodFlag + 1 ));
    double *xold = CreateArray<double>(sysDim);

    //	Set the initial condition.

    for(n=0; n<sysDim; ++n)
	{
        x[n] = solution[n][0];
	}

    //	Forward solve loop.

    for(n=1; n<=nsteps; ++n)
    {

		//	Initialize the current time step variables.

        double tStart = tNodes[n-1];	// Starting time for step.
        double tEnd = tNodes[n];		// Ending time for step.
        double dt = tEnd - tStart;		// Compute the delta t.
		
        //	Store previous step values.

        for (i=0; i<sysDim; ++i)
		{
            xold[i] = solution[i][n-1];
		}

        // We start with the working step equal to the full step and attempt
        // a quasiNewton solve.  If it works, great!  If not, we do a bisection
        // on the step and try it again.

        for ( ;; )
        {

            //  Make a prediction for f at current time plus delta time.
            //	For the higher order method, we make the prediction for the
            //	intermediate node along with the end node.
            //
            //	Last argument = 0 means use previous values.
            //	Last argument = 1

            int preopt = 0;
            predict( xold, xp, tStart, tEnd, sysDim, preopt );
            if ( MethodFlag == Method_DG1 )
            {
                double const tMid = tStart + dt / 3.0;
                predict( xold, &xp[ sysDim ], tStart, tMid, sysDim, preopt );
            }

            //	Do newtons method.
            bool const converged = newton( x, xp, xout, xold, tEnd, tStart, sysDim, n );

            //	If the residual decreased,
            if ( converged )
            {
                if ( tEnd == tNodes[ n ] )   	// We've completed this time step.
                {
                    break;
                }
                else			// We've completed an increment towards the time step.
                {
                    tStart = tEnd;
                    tEnd = tNodes[ n ];
                    dt = tEnd - tStart;
                    for ( int i = 0;i < sysDim; ++i )
                    {
                        xold[ i ] = xp[ i ];
                    }
                    continue;
                }
            }

            //	If we didn't converge, try half the dt value.
            dt *= 0.5;
            tEnd = tStart + dt;

            if ( dt < 1.0e-14 )
            {
				//_fcloseall();
				std::string msg = "forwardSolve: did not converge";
				throw std::runtime_error ( msg );
			}
		}

		//	Store the results.
		for(i=0; i<sysDim; ++i)
		{
			solution[i][n] = xout[i];
		}

	}
    DeleteArray(x);
    DeleteArray(xout);
    DeleteArray(xp);
    DeleteArray(xold);
	return 0;
}

double **GForwardSolve::getSolution()
{
  assert(isSetup_);
	return solution;
}

double **GForwardSolve::getDerivatives()
{
  assert(isSetup_);
	return derivatives;
}

double *GForwardSolve::getMesh()
{
  assert(isSetup_);
	return tNodes;
}

void GForwardSolve::printSolution()
{
  assert(isSetup_);
	FILE *fp;
	fp = fopen("forward_solution.txt","w");
	if(fp==(FILE *) NULL)
	{
		printf("Error opening forward problem output file.\n");
		return;
	}
	for(int i=0;i<nsteps+1;++i)
	{
		fprintf(fp,"%13.10lf ",tNodes[i]);
		for(int n=0;n<sysDim;++n)
		{
			fprintf (fp,"%15.12lf",solution[n][i]);
			if (n<sysDim-1) fprintf (fp," ");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void GForwardSolve::printSolution(int fileID)
{
  assert(isSetup_);
	FILE *fp;
	char buffer[32];
	//sprintf_s(buffer, 32, "forward_solution_%04d.txt", fileID);
	sprintf(buffer, "forward_solution_%04d.txt", fileID);
	fp = fopen(buffer,"w");
	if(fp==(FILE *) NULL)
	{
		printf("Error opening forward problem output file.\n");
		return;
	}
	for(int i=0;i<nsteps+1;++i)
	{
		fprintf(fp,"%13.10lf ",tNodes[i]);
		for(int n=0;n<sysDim;++n)
		{
			fprintf (fp,"%15.12lf",solution[n][i]);
			if (n<sysDim-1) fprintf (fp," ");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void GForwardSolve::predict ( double *yin, double *yp, double const told, 
	       double const tnew, int const n, int const predictor)
{
  assert(isSetup_);
	int i;
	double *y;

	switch(predictor)
	{
	case 0:		// Use the previous value.
		for( i = 0; i < n; ++i )
		{
			yp[i] = yin[i];
		}
		break;
	case 1:		//
		y = new double[n];
		modelPtr->calcDerivs(yin, y, told);
		for( i = 0; i < n; ++i )
		{
			yp[i] = yin[i] + y[i]*(tnew-told)*0.5;
		}
		delete y;
		break;
	}

	return;

}

void GForwardSolve::updateMesh(double *newMesh, int nNodes)
{
  assert(isSetup_);
	DeleteArray(solution,sysDim);
	DeleteArray(derivatives,sysDim);
	DeleteArray(tNodes);
	tNodes = CreateArray<double>(nNodes);
	for(int i=0;i<nNodes;i++)
	{
		tNodes[i] = newMesh[i];
	}
	nsteps = nNodes-1;
	allocateSolutionSpace();
	// Reset the ICs.
    for (int i = 0; i < icCnt; ++i )
    {
        solution[i][0] = IC[i];		// Initial state
    }
}

void GForwardSolve::updateIC(std::vector<double> InitCond)
{
  assert(isSetup_);
    icCnt = InitCond.size();
	IC = CreateArray<double>(icCnt);
    for (int i = 0; i < icCnt; ++i )
    {
		IC[i] = InitCond[i];
        solution[i][0] = InitCond[i];		// Initial state
    }
}

} // namespace GAASP

//--- end of definitions for GForwardSolve ---
