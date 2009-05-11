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
//	File:	  GAdjointSolve.cpp
//	Class:	  GAdjointSolve
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
#include "GAdjointSolve.h"
#include "Spline.h"

namespace GAASP {

using std::numeric_limits;

// ----------------------------------------------------------------------------
//	member constants
// ----------------------------------------------------------------------------

char const * const GAdjointSolve::version =		// class version number
	"1.0.0.0";	// major.minor.release.build

typedef void (*TFunction) ( double *, double *, double );

// ----------------------------------------------------------------------------
//	 constructors and destructor
// ----------------------------------------------------------------------------

GAdjointSolve::GAdjointSolve() 
{ 
  isSetup_ = false;
}

// 	constructors
GAdjointSolve::GAdjointSolve(shared_ptr<GModelBase> model, TAdjointControl adjCon, double **fs, double *fNodes)
{
  Setup(model,adjCon,fs,fNodes);
}

void GAdjointSolve::Setup(shared_ptr<GModelBase> model, TAdjointControl adjCon, double **fs, double *fNodes)
{
	int i, j;
	MethodFlag = adjCon.methodFlag;
	ErrorFlag = adjCon.errorFlag;
	nsteps = int((adjCon.endTime-adjCon.startTime)/adjCon.stepSize);
	ssteps = (adjCon.methodFlag == Method_DG0 ? nsteps : nsteps*3);
	endTime = adjCon.endTime;
	modelPtr = model;
	sysDim = modelPtr->getDim();
	errComp = adjCon.errcomp;
  isSetup_ = true;

	// Storage for the forward solution

	fsoln = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<nsteps+1; ++i)
	{
		for(j = 0; j < sysDim; ++j)
		{
			if(fs != 0)
				fsoln[j][i] = fs[j][i];
			else
				fsoln[j][i] = 0.0;
		}
	}

	// Generate the mesh for the adjoint solution

	tNodes = CreateArray<double>(nsteps+1);
	for(i=0; i<nsteps+1; i++)
	{
		tNodes[i] = adjCon.startTime + i*adjCon.stepSize;
	}

	// Generate the mesh that may include intermediate nodes.

	generateAdjointMesh();

	// Generate the spline for the forward solution.

	y2 = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<sysDim; ++i)
		spline(tNodes, fsoln[i], nsteps+1, y2[i]);

    spln = CreateArray<double>(sysDim, nsteps+1);	// spline results
    for(i=0; i<nsteps+1; ++i)
    {
		for(j=0; j<sysDim; ++j)
		{	
			double yspln = 0.0;
		    splint(tNodes, fsoln[j], y2[j], nsteps + 1, dNodes[i], &yspln);
//		    spln[j][i] = forwardSolution[j][i]; //yspln;
		}
    }

	// Generate the adjoint mesh.

	Initialize ();
}

GAdjointSolve::GAdjointSolve(shared_ptr<GModelBase> model, TAdjointControl adjCon, std::string const &fsolnFN)
{
  Setup(model,adjCon,fsolnFN);
}

void GAdjointSolve::Setup(shared_ptr<GModelBase> model, TAdjointControl adjCon, std::string const &fsolnFN)
{
  isSetup_ = true;
}

GAdjointSolve::~GAdjointSolve ()
{
	Clear ();
}

// ----------------------------------------------------------------------------
//	public functions
// ----------------------------------------------------------------------------

//	Clear
// 	"clear" data members
void GAdjointSolve::Clear ()
{
}

// ----------------------------------------------------------------------------
//	protected functions
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//	private functions
// ----------------------------------------------------------------------------

//	Initialize
void GAdjointSolve::Initialize ()
{
  assert(isSetup_);
	machEpsilon = numeric_limits<double>::epsilon();
}

//	Copy
// 	copy to this
void GAdjointSolve::Copy (
	GAdjointSolve const & object)
{
	if ( &object )
	{
	}
}

int GAdjointSolve::solve()
{
  assert(isSetup_);
	int i;			// Loop counter
    int s;			// Solution storage size (1 for low order, 3 for high)
    int index;		// Index to bypass spline nodes

    double T;		// Ending simulation time
    double **y = 0;		// Solution for the current time.
    double *ytmp = 0;   // Mesh node only

    // Main time loop
    // For METHOD = 1, we need extra storage for the solution on the
    // intermediate RK nodes.

    s = MethodFlag == Method_DG0 ? 1 : 3;
    ytmp = new double[sysDim];

    y = new double * [sysDim];
    for(i=0; i<sysDim; ++i)
        y[i] = new double[s];

    // Set the initial conditions.

	solution = CreateArray<double>(sysDim,ssteps+1);
    for(i=0; i<sysDim; ++i)
        solution[i][ssteps] = y[i][0] = 0.0;

    // Set initial dual data based on psi.

    switch(ErrorFlag)
    {
        case 0:      // Endpoint error, phi(T) = 1 for chosen component
            solution[errComp-1][ssteps] = y[errComp][0] = 1.0;
            break;
        case 1:      // Integral error, phi(T) = 0
        case 2:		 // Weighted by distance
        case 3:
            solution[errComp-1][ssteps] = y[errComp-1][0] = 0.0;
            break;
    }

    // Let t=T-s and use derivative wrt s.  This allows moving forward in
    // time, but we need to use T-s in the arguments.

    T = tNodes[ nsteps ];

    for ( int i = 1; i <= nsteps; ++i )
    {
		// Initialize the current time step variables.
		int iStart = i - 1;		// Index for starting time of current interval
		int iEnd = i;			// Index for ending time of current interval
		double tStart = tNodes[ iStart ]; 	// Starting time for step.
		double tEnd = tNodes[ iEnd ];	// Ending time for step.

		// Perform the solve based on chosen method.
		for ( int j = 0; j < sysDim; ++j )
			ytmp[ j ] = y[ j ][ 0 ];

	    // y above is not used below - replaced by fdMethod - !!! Tom
		for ( int n = 0; n < sysDim; ++n )
			delete [] y[ n ];
		delete [] y;

		y = fdMethod( ytmp, ytmp, ( T - tEnd ), ( T - tStart ), T, nsteps, sysDim );

        // Store the results.

        for ( int j = 0; j < sysDim; ++j )
        {
            switch ( MethodFlag )
            {
                case Method_DG0:
                    solution[ j ][ ssteps - i ] = y[ j ][ 0 ];
                    break;
                case Method_DG1:
                    index = ssteps - i * 3;
                    solution[ j ][ index ] = y[ j ][ 0 ];
                    solution[ j ][ index + 1 ] = y[ j ][ 2 ];
                    solution[ j ][ index + 2 ] = y[ j ][ 1 ];
                    break;
            }
        }
    }

    for ( int i = 0; i < sysDim; ++i )
        delete [] y[ i ];
    delete [] y;
    delete [] ytmp;
	return 0;
}

double **GAdjointSolve::getSolution()
{
  assert(isSetup_);
	return solution;
}

void GAdjointSolve::printSolution()
{
  assert(isSetup_);
	FILE *fp;
	fp = fopen("adjoint_solution.txt","w");
	if(fp==(FILE *) NULL)
	{
		printf("Error opening forward problem output file.\n");
		return;
	}
	for(int i=0; i<ssteps+1; i+=3)
	{
		fprintf(fp,"%13.10lf ",dNodes[i]);
		for(int n=0;n<sysDim;++n)
		{
			fprintf (fp,"%15.12lf",solution[n][i]);
			if (n<sysDim-1) fprintf (fp," ");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void GAdjointSolve::printSolution(int fileID)
{
  assert(isSetup_);
	FILE *fp;
	char buffer[32];
  //sprintf_s(buffer, 32, "adjoint_solution_%04d.txt", fileID);
  sprintf(buffer, "adjoint_solution_%04d.txt", fileID);
	fp = fopen(buffer,"w");
	if(fp==(FILE *) NULL)
	{
		printf("Error opening forward problem output file.\n");
		return;
	}
	for(int i=0; i<ssteps+1; i+=3)
	{
		fprintf(fp,"%13.10lf ",dNodes[i]);
		for(int n=0;n<sysDim;++n)
		{
			fprintf (fp,"%15.12lf",solution[n][i]);
			if (n<sysDim-1) fprintf (fp," ");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void GAdjointSolve::allocateSolutionSpace()
{
  assert(isSetup_);
	solution = CreateArray<double>(sysDim,nsteps+1);
}

void GAdjointSolve::generateAdjointMesh()
{
  assert(isSetup_);

	int i, n;
	double const C1 = 0.2113248654051871;	// .5 - sqrt(3)/6
    double const C2 = 0.7886751345948130;	// .5 + sqrt(3)/6

    ///////////////////////////////////////////////////////////////////////////
    //	Generate the mesh based on type.
    ///////////////////////////////////////////////////////////////////////////

    switch (MethodFlag)
    {
	case Method_DG0:

	    dNodes = CreateArray<double>(nsteps + 1);
	    for(i=0; i<nsteps; ++i)
			dNodes[i] = tNodes[i];
	    break;

    case Method_DG1:  // Copies input mesh and adds interior nodes. (cG(2) adjoint spline)

        int const ssteps = nsteps * 3;
        dNodes = CreateArray<double>(ssteps + 1 );
		n = 0;
		for(i=0; i<nsteps; ++i)  // !!! max = ssteps?
        {
            dNodes[n] = tNodes[i];
			dNodes[n+1] = tNodes[i] + C1 * (tNodes[i+1] - tNodes[i]);
            dNodes[n+2] = tNodes[i] + C2 * (tNodes[i+1] - tNodes[i]);
            n += 3;
        }
        dNodes[ssteps] = tNodes[nsteps];
        break;
    }
}

//--- end of definitions for GAdjointSolve ---

void GAdjointSolve::updateData(double *newMesh, int nSteps, double **fs)
{
  assert(isSetup_);
	int i, j;

	nsteps = nSteps;
	ssteps = (MethodFlag == Method_DG0 ? nsteps : nsteps*3);

	DeleteArray(tNodes);
	tNodes = CreateArray<double>(nsteps+1);
	for(i=0; i<nsteps+1; i++)
	{
		tNodes[i] = newMesh[i];
	}

	generateAdjointMesh();

	DeleteArray(fsoln,sysDim);
	fsoln = CreateArray<double>(sysDim,nSteps+1);
	for(i=0; i<nsteps+1; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			fsoln[j][i] = fs[j][i];
		}
	}
	y2 = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<sysDim; ++i)
		spline(tNodes, fsoln[i], nsteps+1, y2[i]);

    spln = CreateArray<double>(sysDim, nsteps+1);	// spline results
    for(i=0; i<nsteps+1; ++i)
    {
		for(j=0; j<sysDim; ++j)
		{	
			double yspln = 0.0;
		    splint(tNodes, fsoln[j], y2[j], nsteps + 1, dNodes[i], &yspln);
//		    spln[j][i] = forwardSolution[j][i]; //yspln;
		}
    }
}

} // namespace GAASP

