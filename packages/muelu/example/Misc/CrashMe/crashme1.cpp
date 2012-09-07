// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * crashme1.cpp
 *
 *  Created on: 30.10.2011
 *      Author: tobias
 */

#define CRASHME1

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_SegregationAFilterFactory.hpp"
#include "MueLu_SegregationATransferFactory.hpp"

//
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"


std::vector<GlobalOrdinal> DefineG();
void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);

  // custom parameters
  LO maxLevels = 3;

  GO maxCoarseSize=1; //FIXME clp doesn't like long long int
  std::string aggOrdering = "natural";
  int minPerAgg=9;
  int maxNbrAlreadySelected=0;

  // read in problem
  Epetra_Map emap(10201,0,*Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;

  std::cout << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("Condif2Mat.mat",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("Condif2Rhs.mat",emap,ptrf);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);

  // Epetra_CrsMatrix -> Xpetra::Operator
  RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<CrsOperator> crsOp = Teuchos::rcp(new CrsOperator(exA));
  RCP<Operator> Op = Teuchos::rcp_dynamic_cast<Operator>(crsOp);

  // Epetra_Vector -> Xpetra::Vector
  RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  // Epetra_Map -> Xpetra::Map
  const RCP< const Map> map = Xpetra::toXpetra(emap);

  //////////////////////////////////////////////////

  std::vector<GlobalOrdinal> dataindices = DefineG();// DefineM();

  // split maps on current proc
  std::vector<int> gids1;
  std::vector<int> gids2;

  // loop over local map ids on current proc
  for (int i=0; i<emap.NumMyElements(); ++i) {
    if(emap.MyLID(i)) {
      const int gid = emap.GID(i);
      bool gidinmap1 = false;
      for (int j=0; j<Teuchos::as<int>(dataindices.size()); j++) {
        if (dataindices[j] == gid) {
          gidinmap1 = true;
          break;
        }
      }
      if(gidinmap1==true)
        gids1.push_back(gid);
      else
        gids2.push_back(gid);
    }
  }

  int gcount1 = 0;
  int len1 = gids1.size();
  Xpetra::toEpetra(comm)->SumAll(&len1,&gcount1,1);
  Teuchos::RCP<Epetra_Map> emap1 = Teuchos::rcp(new Epetra_Map(gcount1,len1,&gids1[0],0,*Xpetra::toEpetra(comm)));
  int gcount2 = 0;
  int len2 = gids2.size();
  Xpetra::toEpetra(comm)->SumAll(&len2,&gcount2,1);
  Teuchos::RCP<Epetra_Map> emap2 = Teuchos::rcp(new Epetra_Map(gcount2,len2,&gids2[0],0,*Xpetra::toEpetra(comm)));

  // split dofs (=nodes) into several independent subdomains
  // build map extractor
  const RCP< const Xpetra::Map<int, int> > map1 = Xpetra::toXpetra(*emap1);
  const RCP< const Xpetra::Map<int, int> > map2 = Xpetra::toXpetra(*emap2);

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(map1);
  xmaps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(map,xmaps);

  //////////////////////////////////////////////////



  RCP<Hierarchy> H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize(maxCoarseSize);

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);

  RCP<SegregationAFilterFactory> segAfiltered = rcp(new SegregationAFilterFactory("A", NULL, map_extractor));
  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory(segAfiltered));
  dropFact->SetVerbLevel(MueLu::Extreme);
  //dropFact->SetFixedBlockSize(2);
  //RCP<PreDropFunctionConstVal> predrop = rcp(new PreDropFunctionConstVal(0.00001));
  //dropFact->SetPreDropFunction(predrop);
  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory(dropFact));
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
    *out << "aggregate ordering :                    RANDOM" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
    *out << "aggregate ordering :                    GRAPH" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  // build transfer operators
  RCP<NullspaceFactory> nspFact = rcp(new NullspaceFactory()); // make sure that we can keep nullspace!!!
#ifdef CRASHME1
  // with
  // std::string aggOrdering = "natural";
  // int minPerAgg=3;
  // int maxNbrAlreadySelected=0;
  // we get MueLu::SegregationATransferFactory::Build: nnz of Ptent != 1 ???
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact,nspFact));
  RCP<PgPFactory> Pfact = rcp( new PgPFactory(TentPFact) );
  RCP<RFactory> Rfact  = rcp( new GenericRFactory(Pfact));
#elif defined(CRASHME2)
  RCP<TentativePFactory> Pfact = rcp(new TentativePFactory(UCAggFact,nspFact));
  RCP<RFactory> Rfact  = rcp( new TransPFactory(Pfact));
#endif
  RCP<RAPFactory> Acfact = rcp( new RAPFactory(Pfact, Rfact) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);
#if defined(CRASHME1)
  RCP<SegregationATransferFactory> transFact = rcp(new SegregationATransferFactory(TentPFact));
#elif defined(CRASHME2)
  RCP<SegregationATransferFactory> transFact = rcp(new SegregationATransferFactory(Pfact));
#endif
  Acfact->AddTransferFactory(transFact);

  Finest->Keep("Aggregates",UCAggFact.get());
  Finest->Keep("Nullspace",nspFact.get());

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 5);
  ifpackList.set("relaxation: damping factor", (SC) 0.7); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(Xpetra::UseEpetra,ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  Teuchos::ParameterList status;
  status = H->FullPopulate(*Pfact,*Rfact,*Acfact,*SmooFact,0,maxLevels);

  H->SetCoarsestSolver(*SmooFact,MueLu::PRE);


  *out << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  Finest->print(*out);

  RCP<Level> coarseLevel = H->GetLevel(1);
  coarseLevel->print(*out);




  RCP<MultiVector> xLsg = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    H->Iterate(*xRhs,10,*xLsg);

    //xLsg->describe(*out,Teuchos::VERB_EXTREME);
  }

  // print out aggregation information
  for(LocalOrdinal l=0; l<H->GetNumLevels()-1;l++) {
    RCP<Level> level = H->GetLevel((int)l);
    ExportAggregates(level, UCAggFact.get(),comm);
  }

  return EXIT_SUCCESS;
}

std::vector<GlobalOrdinal> DefineG()
{
  int length = 2751;
  GlobalOrdinal data[] = {
      1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1741,1742,1743,1744,1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,1766,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1939,1940,1941,1942,1943,1944,1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,2039,2040,2041,2042,2043,2044,2045,2046,2047,2048,2049,2050,2051,2052,2053,2054,2055,2056,2057,2058,2059,2060,2061,2062,2063,2064,2065,2066,2067,2068,2069,2070,2071,2072,2073,2074,2075,2076,2077,2138,2139,2140,2141,2142,2143,2144,2145,2146,2147,2148,2149,2150,2151,2152,2153,2154,2155,2156,2157,2158,2159,2160,2161,2162,2163,2164,2165,2166,2167,2168,2169,2170,2171,2172,2173,2174,2175,2176,2177,2178,2179,2180,2238,2239,2240,2241,2242,2243,2244,2245,2246,2247,2248,2249,2250,2251,2252,2253,2254,2255,2256,2257,2258,2259,2260,2261,2262,2263,2264,2265,2266,2267,2268,2269,2270,2271,2272,2273,2274,2275,2276,2277,2278,2279,2280,2281,2282,2283,2338,2339,2340,2341,2342,2343,2344,2345,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2361,2362,2363,2364,2365,2366,2367,2368,2369,2370,2371,2372,2373,2374,2375,2376,2377,2378,2379,2380,2381,2382,2383,2384,2385,2386,2438,2439,2440,2441,2442,2443,2444,2445,2446,2447,2448,2449,2450,2451,2452,2453,2454,2455,2456,2457,2458,2459,2460,2461,2462,2463,2464,2465,2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,2476,2477,2478,2479,2480,2481,2482,2483,2484,2485,2486,2487,2488,2489,2539,2540,2541,2542,2543,2544,2545,2546,2547,2548,2549,2550,2551,2552,2553,2554,2555,2556,2557,2558,2559,2560,2561,2562,2563,2564,2565,2566,2567,2568,2569,2570,2571,2572,2573,2574,2575,2576,2577,2578,2579,2580,2581,2582,2583,2584,2585,2586,2587,2588,2589,2590,2591,2639,2640,2641,2642,2643,2644,2645,2646,2647,2648,2649,2650,2651,2652,2653,2654,2655,2656,2657,2658,2659,2660,2661,2662,2663,2664,2665,2666,2667,2668,2669,2670,2671,2672,2673,2674,2675,2676,2677,2678,2679,2680,2681,2682,2683,2684,2685,2686,2687,2688,2689,2690,2691,2692,2693,2739,2740,2741,2742,2743,2744,2745,2746,2747,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2760,2761,2762,2763,2764,2765,2766,2767,2768,2769,2770,2771,2772,2773,2774,2775,2776,2777,2778,2779,2780,2781,2782,2783,2784,2785,2786,2787,2788,2789,2790,2791,2792,2793,2794,2795,2796,2840,2841,2842,2843,2844,2845,2846,2847,2848,2849,2850,2851,2852,2853,2854,2855,2856,2857,2858,2859,2860,2861,2862,2863,2864,2865,2866,2867,2868,2869,2870,2871,2872,2873,2874,2875,2876,2877,2878,2879,2880,2881,2882,2883,2884,2885,2886,2887,2888,2889,2890,2891,2892,2893,2894,2895,2896,2897,2898,2940,2941,2942,2943,2944,2945,2946,2947,2948,2949,2950,2951,2952,2953,2954,2955,2956,2957,2958,2959,2960,2961,2962,2963,2964,2965,2966,2967,2968,2969,2970,2971,2972,2973,2974,2975,2976,2977,2978,2979,2980,2981,2982,2983,2984,2985,2986,2987,2988,2989,2990,2991,2992,2993,2994,2995,2996,2997,2998,2999,3000,3041,3042,3043,3044,3045,3046,3047,3048,3049,3050,3051,3052,3053,3054,3055,3056,3057,3058,3059,3060,3061,3062,3063,3064,3065,3066,3067,3068,3069,3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,3080,3081,3082,3083,3084,3085,3086,3087,3088,3089,3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,3100,3101,3102,3142,3143,3144,3145,3146,3147,3148,3149,3150,3151,3152,3153,3154,3155,3156,3157,3158,3159,3160,3161,3162,3163,3164,3165,3166,3167,3168,3169,3170,3171,3172,3173,3174,3175,3176,3177,3178,3179,3180,3181,3182,3183,3184,3185,3186,3187,3188,3189,3190,3191,3192,3193,3194,3195,3196,3197,3198,3199,3200,3201,3202,3203,3204,3242,3243,3244,3245,3246,3247,3248,3249,3250,3251,3252,3253,3254,3255,3256,3257,3258,3259,3260,3261,3278,3279,3280,3281,3282,3283,3284,3285,3286,3287,3288,3289,3290,3291,3292,3293,3294,3295,3296,3297,3298,3299,3300,3301,3302,3303,3304,3305,3306,3343,3344,3345,3346,3347,3348,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3383,3384,3385,3386,3387,3388,3389,3390,3391,3392,3393,3394,3395,3396,3397,3398,3399,3400,3401,3402,3403,3404,3405,3406,3407,3408,3444,3445,3446,3447,3448,3449,3450,3451,3452,3453,3454,3455,3456,3457,3458,3459,3487,3488,3489,3490,3491,3492,3493,3494,3495,3496,3497,3498,3499,3500,3501,3502,3503,3504,3505,3506,3507,3508,3509,3510,3544,3545,3546,3547,3548,3549,3550,3551,3552,3553,3554,3555,3556,3557,3558,3559,3590,3591,3592,3593,3594,3595,3596,3597,3598,3599,3600,3601,3602,3603,3604,3605,3606,3607,3608,3609,3610,3611,3645,3646,3647,3648,3649,3650,3651,3652,3653,3654,3655,3656,3657,3658,3659,3693,3694,3695,3696,3697,3698,3699,3700,3701,3702,3703,3704,3705,3706,3707,3708,3709,3710,3711,3712,3713,3746,3747,3748,3749,3750,3751,3752,3753,3754,3755,3756,3757,3758,3759,3796,3797,3798,3799,3800,3801,3802,3803,3804,3805,3806,3807,3808,3809,3810,3811,3812,3813,3814,3815,3847,3848,3849,3850,3851,3852,3853,3854,3855,3856,3857,3858,3859,3860,3898,3899,3900,3901,3902,3903,3904,3905,3906,3907,3908,3909,3910,3911,3912,3913,3914,3915,3916,3948,3949,3950,3951,3952,3953,3954,3955,3956,3957,3958,3959,3960,4001,4002,4003,4004,4005,4006,4007,4008,4009,4010,4011,4012,4013,4014,4015,4016,4017,4018,4049,4050,4051,4052,4053,4054,4055,4056,4057,4058,4059,4060,4061,4103,4104,4105,4106,4107,4108,4109,4110,4111,4112,4113,4114,4115,4116,4117,4118,4119,4150,4151,4152,4153,4154,4155,4156,4157,4158,4159,4160,4161,4205,4206,4207,4208,4209,4210,4211,4212,4213,4214,4215,4216,4217,4218,4219,4220,4221,4251,4252,4253,4254,4255,4256,4257,4258,4259,4260,4261,4262,4307,4308,4309,4310,4311,4312,4313,4314,4315,4316,4317,4318,4319,4320,4321,4322,4352,4353,4354,4355,4356,4357,4358,4359,4360,4361,4362,4363,4409,4410,4411,4412,4413,4414,4415,4416,4417,4418,4419,4420,4421,4422,4423,4424,4453,4454,4455,4456,4457,4458,4459,4460,4461,4462,4463,4464,4511,4512,4513,4514,4515,4516,4517,4518,4519,4520,4521,4522,4523,4524,4525,4554,4555,4556,4557,4558,4559,4560,4561,4562,4563,4564,4565,4613,4614,4615,4616,4617,4618,4619,4620,4621,4622,4623,4624,4625,4626,4655,4656,4657,4658,4659,4660,4661,4662,4663,4664,4665,4666,4714,4715,4716,4717,4718,4719,4720,4721,4722,4723,4724,4725,4726,4727,4728,4756,4757,4758,4759,4760,4761,4762,4763,4764,4765,4766,4767,4816,4817,4818,4819,4820,4821,4822,4823,4824,4825,4826,4827,4828,4829,4857,4858,4859,4860,4861,4862,4863,4864,4865,4866,4867,4868,4918,4919,4920,4921,4922,4923,4924,4925,4926,4927,4928,4929,4930,4958,4959,4960,4961,4962,4963,4964,4965,4966,4967,4968,4969,5019,5020,5021,5022,5023,5024,5025,5026,5027,5028,5029,5030,5031,5032,5059,5060,5061,5062,5063,5064,5065,5066,5067,5068,5069,5070,5086,5087,5088,5089,5090,5091,5120,5121,5122,5123,5124,5125,5126,5127,5128,5129,5130,5131,5132,5133,5160,5161,5162,5163,5164,5165,5166,5167,5168,5169,5170,5171,5186,5187,5188,5189,5190,5191,5192,5193,5194,5222,5223,5224,5225,5226,5227,5228,5229,5230,5231,5232,5233,5234,5261,5262,5263,5264,5265,5266,5267,5268,5269,5270,5271,5272,5286,5287,5288,5289,5290,5291,5292,5293,5294,5295,5296,5323,5324,5325,5326,5327,5328,5329,5330,5331,5332,5333,5334,5335,5362,5363,5364,5365,5366,5367,5368,5369,5370,5371,5372,5373,5387,5388,5389,5390,5391,5392,5393,5394,5395,5396,5397,5424,5425,5426,5427,5428,5429,5430,5431,5432,5433,5434,5435,5436,5463,5464,5465,5466,5467,5468,5469,5470,5471,5472,5473,5474,5488,5489,5490,5491,5492,5493,5494,5495,5496,5497,5498,5499,5525,5526,5527,5528,5529,5530,5531,5532,5533,5534,5535,5536,5537,5565,5566,5567,5568,5569,5570,5571,5572,5573,5574,5575,5588,5589,5590,5591,5592,5593,5594,5595,5596,5597,5598,5599,5600,5626,5627,5628,5629,5630,5631,5632,5633,5634,5635,5636,5637,5638,5666,5667,5668,5669,5670,5671,5672,5673,5674,5675,5676,5689,5690,5691,5692,5693,5694,5695,5696,5697,5698,5699,5700,5701,5727,5728,5729,5730,5731,5732,5733,5734,5735,5736,5737,5738,5739,5767,5768,5769,5770,5771,5772,5773,5774,5775,5776,5777,5778,5790,5791,5792,5793,5794,5795,5796,5797,5798,5799,5800,5801,5802,5828,5829,5830,5831,5832,5833,5834,5835,5836,5837,5838,5839,5840,5868,5869,5870,5871,5872,5873,5874,5875,5876,5877,5878,5879,5891,5892,5893,5894,5895,5896,5897,5898,5899,5900,5901,5902,5903,5904,5929,5930,5931,5932,5933,5934,5935,5936,5937,5938,5939,5940,5941,5969,5970,5971,5972,5973,5974,5975,5976,5977,5978,5979,5980,5992,5993,5994,5995,5996,5997,5998,5999,6000,6001,6002,6003,6004,6005,6030,6031,6032,6033,6034,6035,6036,6037,6038,6039,6040,6041,6042,6071,6072,6073,6074,6075,6076,6077,6078,6079,6080,6081,6082,6083,6093,6094,6095,6096,6097,6098,6099,6100,6101,6102,6103,6104,6105,6106,6131,6132,6133,6134,6135,6136,6137,6138,6139,6140,6141,6142,6143,6172,6173,6174,6175,6176,6177,6178,6179,6180,6181,6182,6183,6184,6185,6186,6187,6188,6189,6190,6194,6195,6196,6197,6198,6199,6200,6201,6202,6203,6204,6205,6206,6207,6232,6233,6234,6235,6236,6237,6238,6239,6240,6241,6242,6243,6244,6273,6274,6275,6276,6277,6278,6279,6280,6281,6282,6283,6284,6285,6286,6287,6288,6289,6290,6291,6292,6293,6294,6295,6296,6297,6298,6299,6300,6301,6302,6303,6304,6305,6306,6307,6308,6333,6334,6335,6336,6337,6338,6339,6340,6341,6342,6343,6344,6345,6374,6375,6376,6377,6378,6379,6380,6381,6382,6383,6384,6385,6386,6387,6388,6389,6390,6391,6392,6393,6394,6395,6396,6397,6398,6399,6400,6401,6402,6403,6404,6405,6406,6407,6408,6409,6434,6435,6436,6437,6438,6439,6440,6441,6442,6443,6444,6445,6446,6476,6477,6478,6479,6480,6481,6482,6483,6484,6485,6486,6487,6488,6489,6490,6491,6492,6493,6494,6495,6496,6497,6498,6499,6500,6501,6502,6503,6504,6505,6506,6507,6508,6509,6510,6534,6535,6536,6537,6538,6539,6540,6541,6542,6543,6544,6545,6546,6547,6577,6578,6579,6580,6581,6582,6583,6584,6585,6586,6587,6588,6589,6590,6591,6592,6593,6594,6595,6596,6597,6598,6599,6600,6601,6602,6603,6604,6605,6606,6607,6608,6609,6610,6611,6635,6636,6637,6638,6639,6640,6641,6642,6643,6644,6645,6646,6647,6648,6678,6679,6680,6681,6682,6683,6684,6685,6686,6687,6688,6689,6690,6691,6692,6693,6694,6695,6696,6697,6698,6699,6700,6701,6702,6703,6704,6705,6706,6707,6708,6709,6710,6711,6712,6735,6736,6737,6738,6739,6740,6741,6742,6743,6744,6745,6746,6747,6748,6749,6780,6781,6782,6783,6784,6785,6786,6787,6788,6789,6790,6791,6792,6793,6794,6795,6796,6797,6798,6799,6800,6801,6802,6803,6804,6805,6806,6807,6808,6809,6810,6811,6812,6813,6835,6836,6837,6838,6839,6840,6841,6842,6843,6844,6845,6846,6847,6848,6849,6881,6882,6883,6884,6885,6886,6887,6888,6889,6890,6891,6892,6893,6894,6895,6896,6897,6898,6899,6900,6901,6902,6903,6904,6905,6906,6907,6908,6909,6910,6911,6912,6913,6914,6935,6936,6937,6938,6939,6940,6941,6942,6943,6944,6945,6946,6947,6948,6949,6950,6983,6984,6985,6986,6987,6988,6989,6990,6991,6992,6993,6994,6995,6996,6997,6998,6999,7000,7001,7002,7003,7004,7005,7006,7007,7008,7009,7010,7011,7012,7013,7014,7015,7034,7035,7036,7037,7038,7039,7040,7041,7042,7043,7044,7045,7046,7047,7048,7049,7050,7051,7084,7085,7086,7087,7088,7089,7090,7091,7092,7093,7094,7095,7096,7097,7098,7099,7100,7101,7102,7103,7104,7105,7106,7107,7108,7109,7110,7111,7112,7113,7114,7115,7116,7128,7129,7130,7131,7132,7133,7134,7135,7136,7137,7138,7139,7140,7141,7142,7143,7144,7145,7146,7147,7148,7149,7150,7151,7186,7187,7188,7189,7190,7191,7192,7193,7194,7195,7196,7197,7198,7199,7200,7201,7202,7203,7204,7205,7206,7207,7208,7209,7210,7211,7212,7213,7214,7215,7216,7217,7227,7228,7229,7230,7231,7232,7233,7234,7235,7236,7237,7238,7239,7240,7241,7242,7243,7244,7245,7246,7247,7248,7249,7250,7251,7252,7290,7291,7292,7293,7294,7295,7296,7297,7298,7299,7300,7301,7302,7303,7304,7305,7306,7307,7308,7309,7310,7311,7312,7313,7314,7315,7316,7317,7318,7328,7329,7330,7331,7332,7333,7334,7335,7336,7337,7338,7339,7340,7341,7342,7343,7344,7345,7346,7347,7348,7349,7350,7351,7352,7396,7397,7398,7399,7400,7401,7402,7403,7404,7405,7406,7407,7408,7409,7410,7411,7412,7413,7414,7415,7416,7417,7418,7419,7429,7430,7431,7432,7433,7434,7435,7436,7437,7438,7439,7440,7441,7442,7443,7444,7445,7446,7447,7448,7449,7450,7451,7452,7453,7503,7504,7505,7506,7507,7508,7509,7510,7511,7512,7513,7514,7515,7516,7517,7518,7519,7520,7530,7531,7532,7533,7534,7535,7536,7537,7538,7539,7540,7541,7542,7543,7544,7545,7546,7547,7548,7549,7550,7551,7552,7553,7608,7609,7610,7611,7612,7613,7614,7615,7616,7617,7618,7619,7620,7621,7630,7631,7632,7633,7634,7635,7636,7637,7638,7639,7640,7641,7642,7643,7644,7645,7646,7647,7648,7649,7650,7651,7652,7653,7654,7655,7710,7711,7712,7713,7714,7715,7716,7717,7718,7719,7720,7721,7722,7732,7733,7734,7735,7736,7737,7738,7739,7740,7741,7742,7743,7744,7745,7746,7747,7748,7749,7750,7751,7752,7753,7754,7755,7756,7757,7811,7812,7813,7814,7815,7816,7817,7818,7819,7820,7821,7822,7823,7833,7834,7835,7836,7837,7838,7839,7840,7841,7842,7843,7844,7845,7846,7847,7848,7849,7850,7851,7852,7853,7854,7855,7856,7857,7858,7859,7912,7913,7914,7915,7916,7917,7918,7919,7920,7921,7922,7923,7934,7935,7936,7937,7938,7939,7940,7941,7942,7943,7944,7945,7946,7947,7948,7949,7950,7951,7952,7953,7954,7955,7956,7957,7958,7959,7960,7961,8013,8014,8015,8016,8017,8018,8019,8020,8021,8022,8023,8024,8036,8037,8038,8039,8040,8041,8042,8043,8044,8045,8046,8047,8048,8049,8050,8051,8052,8053,8054,8055,8056,8057,8058,8059,8060,8061,8062,8115,8116,8117,8118,8119,8120,8121,8122,8123,8124,8125,8139,8140,8141,8142,8143,8144,8145,8146,8147,8148,8149,8150,8151,8152,8153,8154,8155,8156,8157,8158,8159,8160,8161,8162,8163,8216,8217,8218,8219,8220,8221,8222,8223,8224,8225,8226,8243,8244,8245,8246,8247,8248,8249,8250,8251,8252,8253,8254,8255,8256,8257,8258,8259,8260,8261,8262,8263,8264,8318,8319,8320,8321,8322,8323,8324,8325,8326,8349,8350,8351,8352,8353,8354,8355,8356,8357,8358,8359,8360,8361,8362,8363,8364,8365,8421,8422,8423,8424,8425,8426,8455,8456,8457,8458,8459,8460,8461,8462,8463,8464,8465,8466,8560,8561,8562,8563,8564,8565,8566
     };
  std::vector<GlobalOrdinal> ret(length,0);
  for (int i=0; i<length; i++)
    ret[i] = data[i];
  return ret;
}

void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm ) {
  if(level->IsAvailable("Aggregates",aggFact) == false) {
    std::cout << "no Aggregates available." << std::endl;
    return;
  }
  if(level->IsAvailable("A") == false || level->IsKept("A")==false) {
    std::cout << "no matrix A available" << std::endl;
    return;
  }
  Teuchos::RCP<Aggregates> aggregates = level->Get< Teuchos::RCP<Aggregates> >("Aggregates",aggFact);
  Teuchos::RCP<Operator> Op = level->Get<Teuchos::RCP<Operator> >("A");

  Teuchos::RCP<LOVector>vertex2AggId_vector = aggregates->GetVertex2AggId();
  Teuchos::RCP<LOVector>procWinner_vector = aggregates->GetProcWinner();
  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

  // 1.) prepare for calculating global aggregate ids
  std::vector<GlobalOrdinal> numAggsGlobal(comm->getSize(),0);
  std::vector<GlobalOrdinal> numAggsLocal(comm->getSize(),0);
  std::vector<GlobalOrdinal> minGlobalAggId(comm->getSize(),0);
  numAggsLocal[comm->getRank()] = aggregates->GetNumAggregates();
  Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, comm->getSize(),&numAggsLocal[0], &numAggsGlobal[0]);
  for(int i=1; i<Teuchos::as<int>(numAggsGlobal.size()); ++i) {
    numAggsGlobal[i] += numAggsGlobal[i-1];
    minGlobalAggId[i] = numAggsGlobal[i-1];
  }

  // 2.) build overlapping vector (global row id (column map) -> global AggId)
  Teuchos::RCP<MultiVector> OverlappingRowId2GlobalAggId = MultiVectorFactory::Build(Op->getColMap(),1);
  Teuchos::ArrayRCP<Scalar> OverlappingRowId2GlobalAggId_data = OverlappingRowId2GlobalAggId->getDataNonConst(0);

  // loop over local aggids
  for(LocalOrdinal lrow=0; lrow<vertex2AggId.size(); ++lrow)
  {
    LocalOrdinal localAggId = vertex2AggId[lrow];
    LocalOrdinal localAggIdProc = procWinner[lrow];
    GlobalOrdinal globalAggId = minGlobalAggId[localAggIdProc] + localAggId;
    OverlappingRowId2GlobalAggId_data[lrow] = globalAggId;
  }

  // 3.) communicate vector (global row id -> global agg id)
  Teuchos::RCP<const Import> importer = ImportFactory::Build(Op->getColMap(), Op->getRowMap());

  Teuchos::RCP<MultiVector> RowMap2GlobalAggId = MultiVectorFactory::Build(Op->getRowMap(),1);
  RowMap2GlobalAggId->doImport(*OverlappingRowId2GlobalAggId,*importer,Xpetra::INSERT);
  Teuchos::ArrayRCP<Scalar> RowMap2GlobalAggId_data = RowMap2GlobalAggId->getDataNonConst(0);

  // 4.) write to file
  std::stringstream filename;
  filename << "aggs_level" << level->GetLevelID() << "_proc" << comm->getRank() << ".out";
  std::ofstream fout(filename.str().c_str());
  for(size_t i=0; i<RowMap2GlobalAggId->getLocalLength(); i++) {
    fout << i << " " << Op->getRowMap()->getGlobalElement(i) << " " << RowMap2GlobalAggId_data[i] << std::endl;
  }
  fout.close();
}

