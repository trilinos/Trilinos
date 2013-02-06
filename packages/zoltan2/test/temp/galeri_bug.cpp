
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraParameters.hpp>

#include <Kokkos_DefaultNode.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using std::string;


template <typename lno_t, typename gno_t, typename scalar_t>
int buildCrsMatrix(int xdim, int ydim, int zdim, string problemType,
                   RCP<const Teuchos::Comm<int> > &comm)
{
  if (comm->getRank() == 0){
    cout << "Create matrix with " << problemType;
    cout << " (and " << xdim;
    if (zdim > 0)
      cout << " x " << ydim << " x " << zdim << " ";
    else if (ydim > 0)
      cout << " x"  << ydim << " x 1 ";
    else
      cout << "x 1 x 1 ";
    cout << " mesh)" << endl;

    cout << "Template Type Sizes:  " << endl
         << "  sizeof(lno_t) = " << sizeof(lno_t) << " " << endl
         << "  sizeof(gno_t) = " << sizeof(gno_t) << " " << endl
         << "  sizeof(scalar_t) = " << sizeof(scalar_t)  << endl;
  }

  Teuchos::CommandLineProcessor tclp;
  Galeri::Xpetra::Parameters<gno_t> params(tclp, xdim, ydim, zdim, problemType);

  if (comm->getRank() == 0) cout << "BUILD MAP" << endl;
  RCP<const Tpetra::Map<lno_t, gno_t> > map =
    rcp(new Tpetra::Map<lno_t, gno_t>(
      params.GetNumGlobalElements(), 0, comm));

  if (comm->getRank() == 0) cout << "BUILD MATRIX" << endl;

  RCP<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> > M_;
  try{
    RCP<Galeri::Xpetra::Problem<Tpetra::Map<lno_t, gno_t>,
                                Tpetra::CrsMatrix<scalar_t, lno_t, gno_t>,
                                Tpetra::MultiVector<scalar_t, lno_t, gno_t> > > Pr=
        Galeri::Xpetra::BuildProblem<scalar_t, lno_t, gno_t,
                                     Tpetra::Map<lno_t, gno_t>,
                                     Tpetra::CrsMatrix<scalar_t, lno_t, gno_t>,
                                     Tpetra::MultiVector<scalar_t, lno_t, gno_t> >
                        (params.GetMatrixType(), map, params.GetParameterList());
    if (comm->getRank() == 0) 
      cout << "AFTER GALERI BuildProblem M_=" << M_ << endl;

    M_ = Pr->BuildMatrix();
    if (comm->getRank() == 0) 
      cout << "AFTER GALERI BuildMatrix  M_=" << M_ << endl;
  }
  catch (exception &e) {    // Probably not enough memory
    cout << "Error returned from Galeri " << e.what() << endl;
    exit(-1);
  }
  if (M_.is_null())
    return 1;
  else 
    return 0;
}

int main(int narg, char **arg)
{
  int ierr, jerr;

  Teuchos::GlobalMPISession mpiSession(&narg, &arg, NULL);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  if (comm->getRank() == 0) cout << "TESTING WITH scalar_t == DOUBLE" << endl;
  ierr = buildCrsMatrix<int, long, double>(10, 10, 10, string("Laplace3D"), comm);
  if (comm->getRank() == 0 && !ierr)
    cout << "SUCCESS for scalar_t == DOUBLE" << endl << endl << endl;

  if (comm->getRank() == 0) cout << "TESTING WITH scalar_t == FLOAT" << endl;
  jerr = buildCrsMatrix<int, long, float>(10, 10, 10, string("Laplace3D"), comm);
  if (comm->getRank() == 0 && !jerr)
    cout << "SUCCESS for scalar_t == FLOAT" << endl << endl << endl;

  if (ierr)
    cout << "FAIL:  M_ was NULL for scalar_t == DOUBLE  " << endl;
  if (jerr)
    cout << "FAIL:  M_ was NULL for scalar_t == FLOAT  " << endl;
  if (!ierr && !jerr)
    cout << "PASS" << endl;
}
