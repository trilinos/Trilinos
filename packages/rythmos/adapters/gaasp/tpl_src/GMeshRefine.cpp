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
//	File:	  GMeshRefine.cpp
//	Class:	  GMeshRefine
//
//	Description:
//	Handles mesh refinement.
// ----------------------------------------------------------------------------
//	Author:	Jeff Sandelin, sandelin@math.colostate.edu, November, 2007
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------

#include <limits>
#include <math.h>
#include "GMeshRefine.h"
#include "Random.h"

namespace GAASP {

using std::numeric_limits;

// ----------------------------------------------------------------------------
//	member constants
// ----------------------------------------------------------------------------

char const * const GMeshRefine::version =		// class version number
	"1.0.0.0";	// major.minor.release.build

// ----------------------------------------------------------------------------
//	 constructors and destructor
// ----------------------------------------------------------------------------
GMeshRefine::GMeshRefine ()
{
  isSetup_ = false;
}

//---- 
GMeshRefine::GMeshRefine (double *oldMesh, double **iErr, TRefinementFlag rFlag,
						  TMethodFlag mFlag, double uTOL, int dim, int nodes)
{
  Setup(oldMesh,iErr,rFlag,mFlag,uTOL,dim,nodes);
}
void GMeshRefine::Setup (double *oldMesh, double **iErr, TRefinementFlag rFlag,
						  TMethodFlag mFlag, double uTOL, int dim, int nodes)
{
	int i, j;
	sysDim = dim;
	userTOL = uTOL;
	nOldNodes = nodes;
	nOldIntervals = nOldNodes - 1;
	oldNodes = new double[nOldNodes];
	intervalError = CreateArray<double>(sysDim,nOldNodes);
	for(i=0; i<nOldNodes; i++)
	{
		oldNodes[i] = oldMesh[i];
	}
	for(i=0; i<nOldNodes; i++)
	{
		for(j=0; j<sysDim; j++)
		{
			if(iErr != 0)
				intervalError[j][i] = iErr[j][i];
			else
				intervalError[j][i] = 0.0;
		}
	}
	rOption = rFlag;
	MethodFlag = mFlag;
  isSetup_ = true;
}

GMeshRefine::~GMeshRefine ()
{
	Clear ();
}

// ----------------------------------------------------------------------------
//	public functions
// ----------------------------------------------------------------------------

//	Clear
// 	"clear" data members
void GMeshRefine::Clear ()
{
}

// ----------------------------------------------------------------------------
//	protected functions
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//	private functions
// ----------------------------------------------------------------------------

//	Initialize
void GMeshRefine::Initialize ()
{
  assert(isSetup_);
	flagsSet = false;
}

//	Copy
// 	copy to this
void GMeshRefine::Copy (
	GMeshRefine const & object)
{
	if ( &object )
	{
	}
}

void GMeshRefine::updateData(double **intEst, double *mesh, int nNodes)
{
  assert(isSetup_);
	int i, j;
	nOldNodes = nNodes;
	nOldIntervals = nOldNodes - 1;
	DeleteArray(oldNodes);
	oldNodes = CreateArray<double>(nNodes);
	intervalError = CreateArray<double>(sysDim,nOldNodes);
	for(i=0; i<nNodes; i++)
	{
		oldNodes[i] = mesh[i];
	}
	for(i=0; i<nNodes; i++)
	{
		for(j=0; j<sysDim; j++)
		{
			intervalError[j][i] = intEst[j][i];
		}
	}
}

void GMeshRefine::refine()
{
  assert(isSetup_);
	computeATOL();
	flagIntervals();
	refineMesh();
}

void GMeshRefine::computeATOL()
{
  assert(isSetup_);
	int i,j;
	double gamma;
	double totalError;
	double totalAbsoluteError;
	double etmp = 0.0;

	totalError = totalAbsoluteError = 0.0;
	for(i=0; i<nOldNodes-1; i++)
	{
		etmp = 0.0;
		for(j=0; j<sysDim; j++)
		{
			etmp += intervalError[j][i];
		}
		totalAbsoluteError += fabs(etmp);
		totalError += etmp;
	}
	gamma = fabs(totalError)/totalAbsoluteError;
	ATOL = userTOL/gamma;
}

void GMeshRefine::flagIntervals()
{
  assert(isSetup_);
	int i, j;				// Loop counters
	int rCount;				// Random number use count
	int rIntervals;			// Number of intervals needing refinement
    int nIntervals;			// Number of negative zone intervals needing refinement
    int pIntervals;			// Number of positive zone intervals needing refinement
    int tIntervals;			// Temporary storage for number of intervals
    int totIntervals;		// Total interval refined - option 5
    int top, bot, mid;		// Nodes in binary search
    int maxRefine;			// Maximum number of intervals to refine
    int maxZoneRefine;		// Maximum number of intervals to refine in current zone
    double newintervals;	// Number of intervals after refinement
    double iTOL;			// Interval tolerance
    double ranInterval;		// Random time
    double ranError;		// Random error level
    double zContrib[2];		// Total contribution for positive and negative zones
    int nzones;				// Number of zones

    int seedFlag = 0;		// Experimental, 0 = fixed seed, 1 = time based
    std::size_t errorSeed = 0;
    std::size_t timeSeed = 0;
	std::vector<double> distTime;
	std::vector<double> distErr;

	flags = CreateArray<int>(nOldNodes-1);

    double fraction;

    int rType = 1;			// 0 = bisection, 1 = equidistribution
    double k, kq;

    int p = 1;							// Default - assume lower order
    if(MethodFlag == Method_DG1) p = 3;	// Higher Order

	int n = nOldNodes-1;				// Set number of intervals
	double const bT = oldNodes[0];		// Set beginning time
    double const eT = oldNodes[n];		// Set ending time

    //	Seed the random number generator object
    if (seedFlag == 0)
    {
		timeSeed = 16384;
		errorSeed = 44730;
    }
    else
    {
		double curtime = time(NULL);
		timeSeed = (std::size_t const) curtime % INT_MAX;
        errorSeed = timeSeed + 387;
    }

    nrel::stats::Random rngTime (timeSeed);
    nrel::stats::Random rngError (errorSeed);

    //	Initialize the refinement flags.
    for(i=0; i<n; ++i)
		flags[i] = 0;

    //	Gather up the interval contributions.
    int *izone = new int[n];					// Positive/negative zone flag
    int *zoneSize = new int[n];					// Number of intervals in zone
    int *zoneType = new int[n];					// Positive (1) or negative (0) zone
    double *zoneContribution = new double[n];	// Total contribution of the zone
    double *intervalContrib = new double[n];	// Interal contribution
    double *cumulError = new double[n];			// Cumulative error
    double totalError = 0.0;					// Sum of interval contributions
    double totAbsError = 0.0;					// Sum of the absolute values of the interval contributions
    double maxErr = 0.0;						// Maximum interval contribution
    for(i=0; i<n; ++i)
    {
		intervalContrib[i] = 0.0;
		for(j=0; j<sysDim; ++j)
		{
			switch (rOption)
			{
			// The first 3 cases all use the same interval error calculation
			case Refine_ByContribution:
			case Refine_Probabilistic1:
			case Refine_Probabilistic2:
				intervalContrib[i] += fabs(intervalError[j][i]);
				break;
			case Refine_WeightedAdaptive:  	//  Need to retain sign.
			    intervalContrib[i] += intervalError[j][i];
			    break;
			}
		}
		totalError += intervalContrib[i];
		totAbsError += fabs(intervalContrib[i]);
		cumulError[i] = fabs(intervalContrib[i]);
		if(i>0)
		    cumulError[i] += cumulError[i-1];
		if(fabs(intervalContrib[i]) > maxErr)
		    maxErr = fabs(intervalContrib[i]);
		izone[i] = intervalContrib[i] < 0 ? 0 : 1;	// Sets zone flags.
    }

    //	Compute the number of refinements.

    newintervals = 0.0;
    for(i=0; i<n; ++i)
		newintervals += pow(fabs(intervalContrib[i] * double(n)/ATOL),(1.0/double(p+1)));
	maxRefine = int(newintervals + 0.5);

    //	Flag intervals for refinement based on refinement option.
    //	Note that options 2, 3, and 4 all refine the same number of intervals,
    //	but the difference is in how the intervals are chosen.

    rIntervals = 0;
    switch(rOption)
    {

	//	Refine all intervals with error contribution > iTOL.
	case Refine_ByContribution:
        iTOL = ATOL / n;
        switch (rType)
        {
        case 0:   	// Bisection
            for (int i=0; i<n; ++i)
            {
                flags[i] = fabs(intervalContrib[i]) < iTOL ? 0 : 1;
                if (flags[i] == 1)
                    rIntervals++;
            }
            break;
        case 1:   	// Equidistribution
            for (int i = 0; i < n; ++i)
            {
                flags[i] = 0;
                if (fabs(intervalContrib[i]) < iTOL) continue;
                kq = ATOL / (n * fabs(intervalContrib[i]));
                kq = pow(kq, double(1.0/p));
                k = (oldNodes[i+1] - oldNodes[i]) * kq;
                flags[i] = int(ceil(1/k));
                rIntervals += flags[i];
            }
            break;
        }
	    break;

	//	Probabilistic refinement based on error contribution.
    case Refine_Probabilistic1:
		// Get random uniform distribution of size maxRefine
		rngError.Uniform(maxRefine, distErr);

		while(rIntervals < maxRefine)
		{
		    ranError = distErr[rIntervals];
		    ranError *= totAbsError;
		    for(i=0; i<n; ++i)
		    {
				if(ranError<cumulError[i] && !flags[i])
				{
				    flags[i] = 1;
				    ++rIntervals;
					break;
				}
		    }
		}
		break;

	//	Probabilistic refinement based on error contribution and time interval.
	case Refine_Probabilistic2:
	
		for(i=0; i<n; ++i)   	// Normalize the errors
		    intervalContrib[i] /= maxErr;

		//	Generate a random time and random error level
		//	Determine the randomly selected interval
		//  Get random uniform distribution of size maxRefine
		rngTime.Uniform(maxRefine, distTime);
		rngError.Uniform(maxRefine, distErr);

		rCount = 0;
		while(rIntervals < maxRefine)
		{
			// Get random numbers
		    ranInterval = distTime[rCount];
		    ranError = distErr[rCount];
			rCount++;

			// If we've used the full set of random numbers, get another set.
			if(rCount == maxRefine)
			{
				rngTime.Uniform(maxRefine, distTime);
				rngError.Uniform(maxRefine, distErr);
				rCount = 0;
			}

		    // Select random interval
		    ranInterval *= (eT - bT);
		    ranInterval += bT;

		    // Select random error level
		    mid = binarySearch(0, n, ranInterval);
		    if(intervalContrib[mid] < ranError) continue;
		    ++flags[mid];
		    ++rIntervals;
		}
		break;

	//	Weighted zone strategy.
    case Refine_WeightedAdaptive:
	
	    //	Compute total contribution from positive and negative zones
	    zContrib[0] = zContrib[1] = 0.0;
	    for(i=0; i<n; ++i)
	        zContrib[izone[i]] += intervalContrib[i];

	    //	Determine number of intervals to refine in + and - zones
	    fraction = fabs(zContrib[0])/(fabs(zContrib[0]) + zContrib[1]);
	    nIntervals = int(maxRefine * fraction + 0.5);
	    pIntervals = int(maxRefine * (1.0-fraction) + 0.5);

	    //	Determine total number of zones, number of intervals in each zone, etc.
	    //	Some variable definitions so I don't have to keep referencing the top of the page.
	    //		nzones			  -  number of zones
	    //		zoneType		  -  0 for negative zone, 1 for positive zone
	    //		zoneSize		  -  number of intervals in zone
	    //		zoneContribution  -  total contribution of zone
	    //		zContrib[i]		  -  total negative (n=0) or positive (n=1) contribution
	    nzones = 1;
	    zoneSize[0] = 1;
	    zoneType[0] = izone[0];
	    zoneContribution[0] = intervalContrib[0];
	    for(i=1; i<n; ++i)
	    {
	        //	Starting a new zone.
	        if(izone[i-1] != izone[i])
	        {
	    		nzones++;
	    		zoneType[nzones-1] = izone[i];
	    		zoneSize[nzones-1] = 1;
	    		zoneContribution[nzones-1] = intervalContrib[i];
	    		continue;
	        }
	        //	Accumulating current zone.
	        zoneSize[nzones-1] ++;
	        zoneContribution[nzones-1] += intervalContrib[i];
	    }
   		
		// Normalize the errors
	    for(i=0; i<n; ++i)
	        intervalContrib[i] /= maxErr;

	    //	Randomly choose intervals within each zone
		std::vector<double> oneErr;
	    bot = totIntervals = 0;
	    for(i=0; i<nzones; ++i)
		{
            //	Determine top end of the zone.
            top = bot + zoneSize[i];

            //	If there's only one interval in the zone, flag it.
            if(zoneSize[i] == 1)
            {
				rngError.Uniform(1, oneErr);
                if(fabs(intervalContrib[bot]) >= oneErr[0])
				{
					flags[bot] = 1;
					totIntervals += 1;
				}
		        bot = top;
                continue;
            }

            //	If there's more than one interval, we choose randomly to determine
            //	how many refinements each interval gets.
            tIntervals = zoneType[ i ] == 0 ? nIntervals : pIntervals;
			maxZoneRefine = int( tIntervals * zoneContribution[ i ] / zContrib[ zoneType[ i ] ] + 0.5 );
            rIntervals = 0;

			// Get random uniform distribution of size maxZoneRefine
			//  - no need to scale since uniform rng object dist is [0,1)
			rngTime.Uniform(maxZoneRefine, distTime);
			rngError.Uniform(maxZoneRefine, distErr);

			rCount = 0;
			while(rIntervals < maxZoneRefine)
            {
				// Get the random numbers to use this round.
				ranInterval = distTime[rCount];
				ranError = distErr[rCount];
				rCount++;

				// If we've used the full set of random numbers, get another set.
				if(rCount == maxZoneRefine)
				{
					rngTime.Uniform(maxZoneRefine, distTime);
					rngError.Uniform(maxZoneRefine, distErr);
					rCount = 0;
				}

				// Select random interval
                ranInterval *= (oldNodes[top] - oldNodes[bot]);
				ranInterval += oldNodes[bot];

				// Select random error level
				mid = binarySearch(bot, top, ranInterval);
                if(fabs(intervalContrib[mid]) < ranError) continue;
				flags[mid] ++;
                rIntervals++;
            }
		    totIntervals += rIntervals;
		    bot = top;
		}
		rIntervals = totIntervals;
		break;

	}	//	End of switch(rOpt)

    delete [] izone;
    delete [] zoneSize;
    delete [] zoneType;
    delete [] zoneContribution;
    delete [] intervalContrib;
    delete [] cumulError;

    nNewIntervals = rIntervals;
}

void GMeshRefine::refineMesh()
{
  assert(isSetup_);
	if(!flagsSet) return;

	int i,j;		//	Loop counters

    int inew;
    int inodes;
    double dt;

    // The method used to flag intervals determines (in part) how intervals are refined.
    // Option 1 uses either bisection or equidistribution of error.
    // Options 3,4,5 have the equidistribution worked into flagIntervals() where
    // the number of intervals to refine is calculated.

    // If we're doing traditional refinement and want equidistribution of error,
    // then we need to do the computations to get the number of new intervals.

    // int rtype = 1;	// 0 = bisection, 1 = equidistribution

    nNewIntervals = nOldIntervals;
    for(i=0; i<nOldIntervals; ++i)
    {
        if(flags[i] == 0)
            continue;
        switch(MethodFlag)
        {
            case Method_DG0:
                nNewIntervals += int(pow(2.0f, flags[i])) - 1;
                break;
            case Method_DG1:
                nNewIntervals += flags[i];
                break;
        }
    }
    newNodes = CreateArray<double>(nNewIntervals + 1);

    //	transfer the old mesh to the new array and insert the
    //	appropriate number of new points.

    inew = 0;
    newNodes[0] = oldNodes[0];
    for(i=1; i<=nOldIntervals; ++i)
    {
        inew++;
        if(flags[i-1] > 0)
        {
			//	Compute the new interval length.
			switch(MethodFlag)
			{
			case Method_DG0:
				inodes = int(pow(2.0f,flags[i-1])) - 1;
				break;
			case Method_DG1:
				inodes = flags[i-1] + 1;
				break;
			}
            dt = (oldNodes[i]-oldNodes[i-1])/inodes;

			for (j=0; j<flags[i-1]; ++j)
            {
                newNodes[inew] = newNodes[inew-1] + dt;
				++inew;
            }
        }
        newNodes[inew] = oldNodes[i];
    }
	nNewNodes = inew + 1;
	return;
}

///////////////////////////////////////////////////////////////////////////////
//	function: binsearch(int bot, int top, double rInt)
//		Searches oldNodes for the required interval.
//		bot	- bottom node
//		top	- top node
//		rInt	- time we want to find
//		mid	- return value
///////////////////////////////////////////////////////////////////////////////
int GMeshRefine::binarySearch(int bot, int top, double const rInt)
{
  assert(isSetup_);
    int ir = -1;	// Search completion flag
    int mid = 0;	// Midpoint

    while(ir<0)
    {
        mid = int((bot + top)/2);
        if (oldNodes[mid] <= rInt && oldNodes[mid+1] > rInt)
            ir = mid;
        else if (rInt < oldNodes[mid])
            top = mid;
        else if (rInt >= oldNodes[mid+1])
            bot = mid;
    }
    return mid;
}

} // namespace GAASP

//--- end of definitions for GMeshRefine ---
