#ifndef INC_nrel_stats_Random_h
#define INC_nrel_stats_Random_h

// ----------------------------------------------------------------------------
//	Organization:	Natural Resource Ecology Laboratory
//			Colorado State University, Fort Collins, CO 80523 USA
//			www.nrel.colostate.edu
//	Project:  Century Soil Organic Matter Model
//	File:	  Random.h
//	Class:	  nrel::stats::Random
//
//	Description:
///	Random number generator using the code from zufall.h/.c.
///	Original description:
///	This package contains a portable random number generator set
///	for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
///	Poisson distributions. The basic module, the uniform generator,
///	uses a lagged Fibonacci series generator:
///
///	              t    = u(n-273) + u(n-607)
///	              u(n) = t - float(int(t))
///
///	where each number generated, u(k), is floating point. Since
///	the numbers are floating point, the left end boundary of the
///	range contains zero. This package is portable except that
///	the test package contains some machine dependent timing data.
///	These are cycle times (in seconds) for NEC SX-3, Fujitsu VP2200,
///	Cray Y-MP, and Sun-4. Select your favorite and comment out the
///	others. There are also vectorization directives for Cray Y-MP
///	machines in the form of "pragma _CRI", which should be ignored
///	although perhaps complained about by other compilers. Otherwise
///	the package is portable and returns the same set of floating
///	point numbers up to word precision on any machine.
///
///	External documentation, "Lagged Fibonacci Random Number Generators
///	for the NEC SX-3," is to be published in the International
///	Journal of High Speed Computing (1994). Otherwise, ask the
///	author:
///		W. P. Petersen
///		IPS, RZ F-5
///		ETHZ
///		CH 8092, Zurich
///		Switzerland
///		e-mail:  wpp@ips.ethz.ch.
// ----------------------------------------------------------------------------
//	Author:	Tom Hilinski, tom.hilinski@colostate.edu, Mar07
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------

#include <vector>

namespace nrel
{
  namespace stats
  {


class Random
{
  public:
	Random (
	  int const useSeed = 1802);		//  0 < seed < 31328

 	// UNIFORM generator functions:
 	// Returns set of n uniforms u[0], ..., u[n-1].
	void Uniform  (int const n, std::vector<double> & u);

 	// NORMAL generator functions:
	// Returns set of n normals g[0], ..., g[n-1] such that
 	// mean <g> = 0, and variance <g**2> = 1.
	void Normal (int const n, std::vector<double> & g);

 	// POISSON generator functions:
	// Returns set of n integers q, with poisson
 	// distribution, density p(q,mu) = exp(-mu) mu**q/q!
	void Poisson (int const n, double const mu, std::vector<int> & q);

  private:

	void InitializeSeedBuffer (int const seed);

 	// Saves seed buffer in saveBuffer[608] for later restarts
 	// using function zufallrs.
	void zufallsv (
	  std::vector<double> & buffer);	//  size = klotz0_1_.size + 1

	// Restores seed buffer in saveBuffer[608] which was previously
	// saved using function zufallsv.
	void zufallrs (
	  std::vector<double> const & buffer);	//  size = klotz0_1_.size + 1

	// Saves seed buffer in saveBuffer[1634] for later restarts
 	// using function normalrs.
	void normalsv (
	  std::vector<double> & buffer);	//  size = 1634

	// Restores seed buffer in saveBuffer[1634] which was previously
 	// saved using function normalsv.
	void normalrs (
	  std::vector<double> const & buffer);	//  size = 1634


	void normal00 ();

	// --------------------------------------------------------------------
	struct klotz0_1_
	{
	    static short const size = 607;
	    std::vector<double> buff;
	    int ptr;

	    klotz0_1_ ()
	      : ptr (0)
	      {
		buff.assign (size, 0.0);
	      }
	};

	klotz0_1_ klotz0_1;

	// --------------------------------------------------------------------
	struct klotz1_1_
	{
	    static short const size = 1024;
	    std::vector<double> xbuff;
	    int first, xptr;

	    klotz1_1_ ()
	      : first (0),
		xptr (0)
	      {
		xbuff.assign (size, 0.0);
	      }
	};

	klotz1_1_ klotz1_1;

};

  } // namespace stats
} // namespace nrel

#endif // INC_nrel_stats_Random_h
