#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_MPISession.hpp"


using namespace Teuchos;

Array<RefCountPtr<Time> > TimeMonitor::timers_;

RefCountPtr<Time> TimeMonitor::getNewTimer(const string& name)
{
	RefCountPtr<Time> rtn = rcp(new Time(name), true);
	timers_.append(rtn);
	return rtn;
}


void TimeMonitor::summarize()
{
	Array<string> names(timers_.length());
	Array<double> timings(timers_.length());

	for (int i=0; i<timers_.length(); i++)
		{
			names[i] = timers_[i]->name();
			timings[i] = timers_[i]->totalElapsedTime();
		}

	int np=1;
	int rank=0;
#ifdef HAVE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (np==1)
		{
			for (int i=0; i<names.length(); i++)
				{
					fprintf(stderr, "%-40s: %g\n", names[i].c_str(), timings[i]);
				}
		}
	else
		{
			Array<double> minTime(timers_.length());
			Array<double> maxTime(timers_.length());
			Array<double> avgTime(timers_.length());
			gatherTimings(timings, minTime, avgTime, maxTime);
			if (rank==0)
				{
					for (int i=0; i<names.length(); i++)
						{
							fprintf(stderr, "%-30s: %-12g %-12g %-12g\n", names[i].c_str(), 
                      minTime[i], avgTime[i], maxTime[i]);
						}
				}
		}

}

void TimeMonitor::gatherTimings(const Array<double>& timings,
                                Array<double>& minTime,
                                Array<double>& avgTime,
                                Array<double>& maxTime)
{
#ifdef HAVE_MPI
	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	void* tPtr = (void*) &(timings[0]);
	void* minPtr = (void*) &(minTime[0]);
	void* avgPtr = (void*) &(avgTime[0]);
	void* maxPtr = (void*) &(maxTime[0]);

	int count = (int) timings.length();

	MPI_Allreduce(tPtr, minPtr, count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(tPtr, avgPtr, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(tPtr, maxPtr, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for (int i=0; i<avgTime.length(); i++)
		{
			avgTime[i] = avgTime[i]/((double) np);
		}
#endif
}

