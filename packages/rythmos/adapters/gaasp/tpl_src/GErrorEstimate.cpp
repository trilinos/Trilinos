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
//	File:	  GErrorEstimate.cpp
//	Class:	  GErrorEstimate
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
#include "GErrorEstimate.h"
#include "Spline.h"
#include "Projections.h"

namespace GAASP {
using std::numeric_limits;

// ----------------------------------------------------------------------------
//	member constants
// ----------------------------------------------------------------------------

char const * const GErrorEstimate::version =		// class version number
	"1.0.0.0";	// major.minor.release.build

typedef void (*TFunction) ( double *, double *, double );

// ----------------------------------------------------------------------------
//	 constructors and destructor
// ----------------------------------------------------------------------------
GErrorEstimate::GErrorEstimate ()
{
  isSetup_ = false;
}

//---- Forward solution, adjoint solution, time nodes, dimension
GErrorEstimate::GErrorEstimate (double **fs, double **as, double **sd, double *nodes, 
								int dim, int steps, shared_ptr<GModelBase> model, TMethodFlag method)
{
  Setup(fs,as,sd,nodes,dim,steps,model,method);
}
void GErrorEstimate::Setup (double **fs, double **as, double **sd, double *nodes, 
								int dim, int steps, shared_ptr<GModelBase> model, TMethodFlag method)
{
	int i, j;
	nsteps = steps;
	MethodFlag = method;
	modelPtr = model;
	sysDim = modelPtr->getDim();
  isSetup_ = true;

	//	Store the forward solution
	fsoln = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<nsteps+1; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			if(fs != 0)
				fsoln[j][i] = fs[j][i];
			else
				fsoln[j][i] = 0.0;
		}
	}

	//	Store the time mesh
	tNodes = CreateArray<double>(nsteps+1);
	for(i=0; i<nsteps+1; i++)
	{
		tNodes[i] = nodes[i];
	}

	//	Generate the adjoint mesh

	generateAdjointMesh();

	//	Store the adjoint solution
	dsoln = CreateArray<double>(sysDim,ssteps+1);
	for(i=0; i<ssteps+1; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			if(as != 0)
				dsoln[j][i] = as[j][i];
			else
				dsoln[j][i] = 0.0;
		}
	}

	//	Store the forward solution derivatives
	sder = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<nsteps+1; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			sder[j][i] = 0.0;
			if(MethodFlag==Method_DG0) continue;
			if(sd == 0) continue;
			sder[j][i] = sd[j][i];
		}
	}


	//	Generate spline of the adjoint solution
	if(as != 0)
	{
		dualy2 = CreateArray<double>(sysDim, ssteps + 1);
		for (i=0; i<sysDim; ++i)
			spline(sNodes, dsoln[i], ssteps+1, dualy2[i]);
	}

	//	Gauss points and weights for 5 point Gauss rule
	t[0] = -0.90617985;
	t[1] = -0.53846931; 
	t[2] = 0.0; 
	t[3] = 0.53846931; 
	t[4] = 0.90617985;
	w[0] = w[4] = 0.23692689;
	w[1] = w[3] = 0.47862867;
	w[2] = 0.56888889;

	Initialize();
}

GErrorEstimate::~GErrorEstimate ()
{
	Clear ();
}

// ----------------------------------------------------------------------------
//	public functions
// ----------------------------------------------------------------------------

//	Clear
// 	"clear" data members
void GErrorEstimate::Clear ()
{
}

// ----------------------------------------------------------------------------
//	protected functions
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//	private functions
// ----------------------------------------------------------------------------

//	Initialize
void GErrorEstimate::Initialize ()
{
  assert(isSetup_);
	allocateSolutionSpace();
}

//	Copy
// 	copy to this
void GErrorEstimate::Copy (
	GErrorEstimate const & object)
{
	if ( &object )
	{
	}
}

void GErrorEstimate::generateAdjointMesh()
{
  assert(isSetup_);
	int i, j;
    double C1 = 0.2113248654051871;	// .5 - sqrt(3)/6
    double C2 = 0.7886751345948130;	// .5 + sqrt(3)/6

	switch(MethodFlag)
	{
	case Method_DG0:
		sNodes = tNodes;
		ssteps = nsteps;
		break;
	case Method_DG1:
		sNodes = CreateArray<double>(nsteps*3+1);
	    j = 0;
		for(i=0; i<nsteps; ++i)
		{
			sNodes[j] = tNodes[i];
			sNodes[j+1] = tNodes[i] + C1 * (tNodes[i+1]-tNodes[i]);
		    sNodes[j+2] = tNodes[i] + C2 * (tNodes[i+1]-tNodes[i]);
			j += 3;
		}
		ssteps = nsteps * 3;
		sNodes[ssteps] = tNodes[nsteps];
		break;
	}
}

void GErrorEstimate::allocateSolutionSpace()
{
  assert(isSetup_);
	intervalError = CreateArray<double>(sysDim,nsteps+1);
	discError = CreateArray<double>(sysDim);
	quadError = CreateArray<double>(sysDim);
}

int GErrorEstimate::compute()
{
  assert(isSetup_);
	int i, j;
	double dt, c, m;
	
	for (i=0; i<nsteps; ++i) 		// Loop through mesh
    {
	    // Compute the slope and intercept for the map from [a,b] to (-1,1)
	    dt = tNodes[i+1] - tNodes[i];
		c = (tNodes[i+1] + tNodes[i])/2.0;
		m = (tNodes[i+1] - tNodes[i])/2.0;
		// Compute the Radau points for this interval
	    tau[0] = tNodes[i] + dt/3.0;
		tau[1] = tNodes[i+1];
		dtau = tau[1] - tau[0];
		// Compute the error contributions
		discretizationError(i,dt,c,m);
		quadratureError(i,dt,c,m);
        for(j=0; j<sysDim; ++j)
		{
			intervalError[j][i] = discError[j] + quadError[j];
		}
    }
	return 0;
}

double GErrorEstimate::getErrorEstimate()
{
  assert(isSetup_);
	int i, j;
	double estimate = 0.0;
    for(i=0; i<nsteps; ++i)
    {
        for(j=0; j<sysDim; ++j)
        {
			estimate += intervalError[j][i];
        }
    }
	return estimate;
}

double **GErrorEstimate::getIntervalContributions()
{
  assert(isSetup_);
	return intervalError;
}

void GErrorEstimate::updateForwardSolution()
{
}

void GErrorEstimate::updateForwardDerivatives()
{
}

void GErrorEstimate::updateAdjointSolution()
{
}

void GErrorEstimate::updateMesh()
{
}

void GErrorEstimate::updateData(double *newMesh, int nNodes, double **fs, double **fd, double **as)
{
  assert(isSetup_);
	int i, j;

	nsteps = nNodes-1;

	DeleteArray(intervalError,sysDim);
	DeleteArray(discError);
	DeleteArray(quadError);
	DeleteArray(tNodes);

	tNodes = CreateArray<double>(nNodes);
	for(i=0; i<nNodes; i++)
	{
		tNodes[i] = newMesh[i];
	}

	DeleteArray(sNodes);
	generateAdjointMesh();

	DeleteArray(fsoln,sysDim);
	fsoln = CreateArray<double>(sysDim,nsteps+1);
	for(i=0; i<nNodes; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			fsoln[j][i] = fs[j][i];
		}
	}
	
	DeleteArray(dsoln,sysDim);
	dsoln = CreateArray<double>(sysDim,ssteps+1);
	for(i=0; i<ssteps+1; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			dsoln[j][i] = as[j][i];
		}
	}

	DeleteArray(dualy2,sysDim);
	dualy2 = CreateArray<double>(sysDim, ssteps + 1);
	for (i=0; i<sysDim; ++i)
		spline(sNodes, dsoln[i], ssteps+1, dualy2[i]);
		
	DeleteArray(sder,sysDim);
	sder = CreateArray<double>(sysDim,nNodes);
	for(i=0; i<nNodes; ++i)
	{
		for(j=0; j<sysDim; ++j)
		{
			sder[j][i] = 0.0;
			if(MethodFlag==Method_DG0) continue;
			sder[j][i] = fd[j][i];
		}
	}

	allocateSolutionSpace();
}

void GErrorEstimate::printErrorEstimates()
{
  assert(isSetup_);
	int i, j;
	FILE *fp;
	fp = fopen("error_contributions.txt","w");
	if(fp==(FILE *) NULL)
	{
		printf("Error opening error estimate output file.\n");
		return;
	}

	fprintf( fp, "\nInterval contributions:\n\n" );

    for(i=0; i<nsteps; ++i)
    {
		fprintf(fp,"%16.14lf",tNodes[i+1]);
        for(j=0; j<sysDim; ++j)
        {
			fprintf(fp,"  %16.14lf",intervalError[j][i]);
        }
        fprintf(fp,"\n");
    }
}

void GErrorEstimate::discretizationError(int iNode, double dt, double c, double m)
{
  assert(isSetup_);
	int i, j;
	int itmp, ijmp;
	double *If, *Iy, *yb, *yout;
	double tt, phi, piphi, uprime;

	If = CreateArray<double>(sysDim);
	Iy = CreateArray<double>(sysDim);
	yb = CreateArray<double>(sysDim);
	yout = CreateArray<double>(sysDim);

	//	Initialize the arrays.
	for(i=0; i<sysDim; ++i)
		Iy[i] = If[i] = 0.0;

	//	Set the node indecies based on method.
	switch(MethodFlag)
	{
	case Method_DG0:
		itmp = iNode;
		ijmp = 0;
		break;
	case Method_DG1:
		itmp = iNode * 3;
		ijmp = 3;
		break;
	}

	//	Evaluate the integrals.
	for(i=0; i<gaussRule; ++i)
	{
        tt = c + m * t[i];
        for(j=0; j<sysDim; ++j)
        {
            splint(sNodes, dsoln[j], dualy2[j], ssteps, tt, &phi);
            piphi = projection(tt, dsoln[j][itmp], dsoln[j][itmp + ijmp],
                                sNodes[itmp], sNodes[itmp + ijmp]);
            uprime = sder[j][iNode + 1];
		    //	INT(yprime,phi-PHI)dt
            Iy[j] += uprime * (phi - piphi) * w[i] * m;
		}
	}
	
	for(i=0; i<gaussRule; ++i)
    {
        tt = c + m * t[i];
        ybar(iNode, tt, yb);
		modelPtr->calcDerivs(yb, yout, tt);
        for(j=0; j<sysDim; ++j)
        {
            splint(sNodes, dsoln[j], dualy2[j], ssteps, tt, &phi);
            piphi = projection(tt, dsoln[j][itmp], dsoln[j][itmp + ijmp],
                                sNodes[itmp], sNodes[itmp + ijmp]);
            uprime = sder[j][iNode + 1];
		    //	INT(f(Y),phi-PHI)dt
            If[j] += yout[j] * (phi - piphi) * w[i] * m;
        }
    }

    //	Compute total discretization error.
    for(i=0; i<sysDim; ++i)
        discError[i] = -Iy[i] + If[i];

	DeleteArray(If);
	DeleteArray(Iy);
	DeleteArray(yb);
	DeleteArray(yout);
}

void GErrorEstimate::quadratureError(int iNode, double dt, double c, double m)
{
  assert(isSetup_);
	int i, j;
	int itmp, ijmp;
	double *If, *Iy, *yb, *yout, *fB;
	double tt, piphi;

	If = CreateArray<double>(sysDim);
	Iy = CreateArray<double>(sysDim);
	yb = CreateArray<double>(sysDim);
	yout = CreateArray<double>(sysDim);
	fB = CreateArray<double>(sysDim);

	//	Initialize the arrays.
	for(i=0; i<sysDim; ++i)
        If[i] = Iy[i] = 0.0;

	//	Set the node indecies based on method.
	switch(MethodFlag)
	{
	case Method_DG0:
		itmp = iNode;
		ijmp = 0;
		break;
	case Method_DG1:
		itmp = iNode * 3;
		ijmp = 3;
		break;
	}

	//	Evaluate the integrals.
    for(i=0; i<gaussRule; ++i)
    {
        tt = c + m * t[i];
        ybar(iNode, tt, yb);
        fBAR(iNode, tt, fB);
		modelPtr->calcDerivs(yb, yout, tt);
        for (j=0; j<sysDim; ++j)
        {
            piphi = projection(tt, dsoln[j][itmp], dsoln[j][itmp + ijmp],
                                sNodes[itmp], sNodes[itmp + ijmp]);
		    //	INT(f(Y),PHI)dt
            If[j] += yout[j] * piphi * w[i] * m;
		    //	INT([f(Y),PHI]_bar)dt
            Iy[ j ] += fB[ j ] * piphi * w[ i ] * m;
        }
    }

    //	Compute total quadrature error.
    for (i=0; i<sysDim; ++i)
        quadError[i] = If[i] - Iy[i];

	DeleteArray(If);
	DeleteArray(Iy);
	DeleteArray(yb);
	DeleteArray(yout);
	DeleteArray(fB);
}

//---- The following functions are utility functions used in the error computations.

//	Evaluation of basis function.
double GErrorEstimate::gamma(int k, double t)
{
  assert(isSetup_);
    double rVal = 0;
    switch(k)
    {
        case 1:
            rVal = (t-tau[1]) / (tau[0]-tau[1]);
            break;
        case 2:
            rVal = (t-tau[0]) / (tau[1]-tau[0]);
            break;
    }
    return rVal;
}

//	Evaluation of Y at the Radau points.
void GErrorEstimate::yrad(int k, int inode, double *rVal)
{
  assert(isSetup_);
    for(int i=0; i<sysDim; ++i)
    {
        switch (k)
        {
            case 0:
                rVal[i] = fsoln[i][inode + 1] - sder[i][inode + 1] * dtau;
                break;
            case 1:
                rVal[i] = fsoln[i][inode + 1];
                break;
        }
    }
}

//	ybar(t) = Y(tau1)*gamma1(t)+Y(tau2)*gamma2(t)
void GErrorEstimate::ybar(int inode, double t, double *rVal)
{
  assert(isSetup_);
    double *ytau1 = new double[sysDim];
    double *ytau2 = new double[sysDim];

    yrad(0, inode, ytau1);
    yrad(1, inode, ytau2);

    for(int i=0; i<sysDim; ++i)
    {
        rVal[i] = ytau1[i] * gamma(1, t) + ytau2[i] * gamma(2, t);
    }

    delete [] ytau1;
    delete [] ytau2;
}

//	fbar(t) = f(Y(tau1),tau1)*gamma1(t)+f(Y(tau2),tau2)*gamma2(t)
void GErrorEstimate::fbar(int inode, double t, double *rVal)
{
  assert(isSetup_);
    double *ytau1 = new double[sysDim];
    double *ytau2 = new double[sysDim];
    double *yout1 = new double[sysDim];
    double *yout2 = new double[sysDim];

    yrad(0,inode, ytau1);
    yrad(1,inode, ytau2);

	modelPtr->calcDerivs(ytau1, yout1, tau[0]);
	modelPtr->calcDerivs(ytau2, yout2, tau[1]);

    for(int i=0; i<sysDim; ++i)
    {
        rVal[i] = yout1[i] * gamma(1, t) + yout2[i] * gamma(2, t);
    }

    delete [] ytau1;
    delete [] ytau2;
    delete [] yout1;
    delete [] yout2;
}

//	fBAR(t) = f(ybar(tau1),tau1)*gamma1(t)+f(ybar(tau2),tau2)*gamma2(t)
void GErrorEstimate::fBAR(int inode, double t, double *rVal)
{
  assert(isSetup_);
    double * yb1 = new double[sysDim];
    double *yb2 = new double[sysDim];
    double *yout1 = new double[sysDim];
    double *yout2 = new double[sysDim];

    ybar(inode, tau[0], yb1);
    ybar(inode, tau[1], yb2);

	modelPtr->calcDerivs(yb1, yout1, tau[0]);
	modelPtr->calcDerivs(yb2, yout2, tau[1]);

    fbar(inode, tau[0], yout1);
    fbar(inode, tau[1], yout2);

    for(int i=0; i<sysDim; ++i)
    {
        rVal[i] = yout1[i] * gamma(1, t) + yout2[i] * gamma(2, t);
    }

    delete [] yb1;
    delete [] yb2;
    delete [] yout1;
    delete [] yout2;
}

} // namespace GAASP

//--- end of definitions for GErrorEstimate ---
