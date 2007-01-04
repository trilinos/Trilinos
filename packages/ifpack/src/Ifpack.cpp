#include "Ifpack_ConfigDefs.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_SPARSKIT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Ifpack_Amesos.h"
#endif
#include "Ifpack_Chebyshev.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StringToIntMap.hpp"

namespace {

const Teuchos::StringToIntMap
precTypeNameToIntMap(
  "parameter \"Prec Type\"", Ifpack::numPrecTypes, Ifpack::precTypeNames
  );

} // namespace

//==============================================================================
const Ifpack::EPrecType Ifpack::precTypeValues[Ifpack::numPrecTypes] =
{
  POINT_RELAXATION
  ,POINT_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION
  ,BLOCK_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION_STAND_ALONE_ILU
#ifdef HAVE_IFPACK_AMESOS
  ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
  ,BLOCK_RELAXATION_AMESOS
  ,AMESOS
  ,AMESOS_STAND_ALONE
#endif // HAVE_IFPACK_AMESOS
  ,IC
  ,IC_STAND_ALONE
  ,ICT
  ,ICT_STAND_ALONE
  ,ILU
  ,ILU_STAND_ALONE
  ,ILUT
  ,ILUT_STAND_ALONE
#ifdef HAVE_IFPACK_SPARSKIT
  ,SPARSKIT
#endif // HAVE_IFPACK_SPARSKIT
  ,CHEBYSHEV
};

//==============================================================================
const char* Ifpack::precTypeNames[Ifpack::numPrecTypes] =
{
  "point relaxation"
  ,"point relaxation stand-alone"
  ,"block relaxation"
  ,"block relaxation stand-alone"
  ,"block relaxation stand-alone (ILU)"
#ifdef HAVE_IFPACK_AMESOS
  ,"block relaxation stand-alone (Amesos)"
  ,"block relaxation (Amesos)"
  ,"Amesos"
  ,"Amesos stand-alone"
#endif
  ,"IC"
  ,"IC stand-alone"
  ,"ICT"
  ,"ICT stand-alone"
  ,"ILU"
  ,"ILU stand-alone"
  ,"ILUT"
  ,"ILUT stand-alone"
#ifdef HAVE_IFPACK_SPARSKIT
  ,"SPARSKIT"
#endif
  ,"Chebyshev"
};

//==============================================================================
const bool Ifpack::supportsUnsymmetric[Ifpack::numPrecTypes] =
{
  true // point relaxation
  ,true // point relaxation stand-alone
  ,true // block relaxation
  ,true // block relaxation stand-alone
  ,true // block relaxation stand-alone (ILU)
#ifdef HAVE_IFPACK_AMESOS
  ,true // block relaxation stand-alone (Amesos)
  ,true // block relaxation (Amesos)
  ,true // Amesos
  ,true // Amesos stand-alone 
#endif
  ,false // IC
  ,false // IC stand-alone
  ,false // ICT
  ,false // ICT stand-alone
  ,true // ILU
  ,true // ILU stand-alone
  ,true // ILUT
  ,true // ILUT stand-alone
#ifdef HAVE_IFPACK_SPARSKIT
  ,true // SPARSKIT
#endif
  ,false // CHEBYSHEV
};

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(EPrecType PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap)
{
  switch(PrecType) {
    case POINT_RELAXATION:
      return(new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(Matrix, Overlap));
    case POINT_RELAXATION_STAND_ALONE:
      return(new Ifpack_PointRelaxation(Matrix));
    case BLOCK_RELAXATION:
      return(new Ifpack_AdditiveSchwarz<
             Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,Overlap));
    case BLOCK_RELAXATION_STAND_ALONE:
      return(new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix));
    case BLOCK_RELAXATION_STAND_ALONE_ILU:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_ILU> >(Matrix));
#ifdef HAVE_IFPACK_AMESOS
    case BLOCK_RELAXATION_STAND_ALONE_AMESOS:
      return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >(Matrix));
    case BLOCK_RELAXATION_AMESOS:
      return(new Ifpack_AdditiveSchwarz<
             Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
    case AMESOS:
      return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
    case AMESOS_STAND_ALONE:
      return(new Ifpack_Amesos(Matrix));
#endif
    case IC:
      return(new Ifpack_AdditiveSchwarz<Ifpack_IC>(Matrix,Overlap));
    case IC_STAND_ALONE:
      return(new Ifpack_IC(Matrix));
    case ICT:
      return(new Ifpack_AdditiveSchwarz<Ifpack_ICT>(Matrix,Overlap));
    case ICT_STAND_ALONE:
      return(new Ifpack_ICT(Matrix));
    case ILU:
      return(new Ifpack_AdditiveSchwarz<Ifpack_ILU>(Matrix,Overlap));
    case ILU_STAND_ALONE:
      return(new Ifpack_ILU(Matrix));
    case ILUT:
      return(new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(Matrix,Overlap));
    case ILUT_STAND_ALONE:
      return(new Ifpack_ILUT(Matrix));
#ifdef HAVE_IFPACK_SPARSKIT
    case SPARSKIT:
      return(new Ifpack_SPARSKIT(Matrix));
#endif
    case CHEBYSHEV:
      return(new Ifpack_Chebyshev(Matrix));
    default:
      TEST_FOR_EXCEPT(true);
      // The only way to get here is if some code developer does a cast like
      // (EPrecType)(anyNumber).  You can never get here by passing in a
      // string value for the preconditioner!
  } // end switch
  return 0; // This will never ever be executed!
}

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(const string PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap)
{
  try {
    return Ifpack::Create(Teuchos::get<EPrecType>(::precTypeNameToIntMap,PrecType),Matrix,Overlap);
  }
  catch( const Teuchos::StringToIntMap::DoesNotExist &excpt ) {
    // The old implementation of this function just silently returned a NULL
    // when a preconditiner type name was not recognized.  If you like this
    // behavior then you should use this function.  If you do not like this
    // behavior, then consider using the Ifpack/Thyra adapter
    // Thyra::IfpackPreconditionerFactory or better yet the Stratimikos
    // wrapper class Thyra::DefaultRealLinearSolverBuilder.
  }
  return 0;
}

// ======================================================================
int Ifpack::SetParameters(int argc, char* argv[],
                          Teuchos::ParameterList& List, string& PrecType,
                          int& Overlap)
{
  // THIS FUNCTION IS VERY INCOMPLETE...

  Teuchos::CommandLineProcessor CLP;

  // prec type
  string ifp_prec_type = "ILU";
  CLP.setOption("ifp-prec-type",&ifp_prec_type,"Preconditioner type");
  // overlap among the processors
  int ifp_overlap = 0;
  CLP.setOption("ifp-overlap",&ifp_overlap,"Overlap among processors");
  // relaxation type
  string ifp_relax_type = "Jacobi";
  CLP.setOption("ifp-relax-type",&ifp_relax_type,"Relaxation type");
  // sweeps (for relax only)
  int ifp_relax_sweeps = 1;
  CLP.setOption("ifp-relax-sweeps",
                &ifp_relax_sweeps,"Number of sweeps for relaxation");
  // damping (for relax only)
  double ifp_relax_damping = 1.0;
  CLP.setOption("ifp-relax-damping",
                &ifp_relax_damping,"Damping for relaxation");
  // partitioner type (for block relaxation only)
  string ifp_part_type = "greedy";
  CLP.setOption("ifp-part-type",&ifp_part_type,"Partitioner type");
  // number of local parts (for block relaxation only)
  int ifp_part_local = 1;
  CLP.setOption("ifp-part-local",&ifp_part_local,"number of local partitions");

  // allow users to specify other options for other packages
  CLP.recogniseAllOptions(false);
  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  // I cannot really set those in the List, I pass them back to the user
  PrecType = ifp_prec_type;
  Overlap = ifp_overlap;

  // set the list here
  List.set("relaxation: type", ifp_relax_type);
  List.set("relaxation: sweeps", ifp_relax_sweeps);
  List.set("relaxation: damping factor", ifp_relax_damping);
  List.set("partitioner: type", ifp_part_type);
  List.set("partitioner: local parts", ifp_part_local);

  return(0);
}
