#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_TEUCHOS)
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#include "Teuchos_CommandLineProcessor.hpp"
#ifdef HAVE_IFPACK_AMESOS
#include "Ifpack_Amesos.h"
#endif

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(const string PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap)
{

  if (PrecType == "point relaxation") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(Matrix, Overlap));
  }
  else if (PrecType == "point relaxation stand-alone") {
    return(new Ifpack_PointRelaxation(Matrix));
  }
  else if (PrecType == "block relaxation") {
    return(new Ifpack_AdditiveSchwarz<
            Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,Overlap));
  }
  else if (PrecType == "block relaxation stand-alone") {
    return(new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(Matrix));
  }
#ifdef HAVE_IFPACK_AMESOS
  else if (PrecType == "block relaxation stand-alone (ILU)") {
    return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >(Matrix));
  }
  else if (PrecType == "block relaxation stand-alone (Amesos)") {
    return(new Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> >(Matrix));
  }
  else if (PrecType == "block relaxation (Amesos)") {
    return(new Ifpack_AdditiveSchwarz<
            Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
  }
  else if (PrecType == "Amesos") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
  }
  else if (PrecType == "Amesos stand-alone" || PrecType == "LU") {
    return(new Ifpack_Amesos(Matrix));
  }
#endif
  else if (PrecType == "IC") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_IC>(Matrix,Overlap));

  } 
  else if (PrecType == "IC stand-alone") {
    return(new Ifpack_IC(Matrix));

  } 
  else if (PrecType == "ICT") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_ICT>(Matrix,Overlap));

  } 
  else if (PrecType == "ICT stand-alone") {
    return(new Ifpack_ICT(Matrix));

  } 
  else if (PrecType == "ILU") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_ILU>(Matrix,Overlap));
  }
  else if (PrecType == "ILU stand-alone") {
    return(new Ifpack_ILU(Matrix));
  }
  else if (PrecType == "ILUT") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(Matrix,Overlap));
  }
  else if (PrecType == "ILUT stand-alone") {
    return(new Ifpack_ILUT(Matrix));
  }
  else
    // nothing understandable
    return(0);

}

// ======================================================================
int Ifpack::SetParameters(int argc, char* argv[],
                          Teuchos::ParameterList& List, string& PrecType,
                          int& Overlap)
{

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
#endif // HAVE_IFPACK_TEUCHOS
