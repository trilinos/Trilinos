#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <iostream>
#include <TpetraExt_MatrixMatrix_def.hpp>

using namespace Teuchos;

typedef Tpetra::Vector<>::local_ordinal_type scalar_type;
typedef Tpetra::Map<> map_type;
typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
typedef Tpetra::Vector<>::mag_type magnitude_type;
typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> crs_matrix_type;
typedef Tpetra::CrsGraph<> crs_graph_type;


template < class TpetraOperatorType>
class PowerMethod{

public:
  typedef typename TpetraOperatorType::scalar_type scalar_type;
  typedef typename TpetraOperatorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraOperatorType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraOperatorType::node_type node_type;

  typedef typename Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vec_type;

  typedef typename vec_type::mag_type mag_type;
  
  static scalar_type run(const TpetraOperatorType & A, const int niters, const mag_type tolerance, std::ostream &out)
  {
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    // A * q = z
    vec_type q(A.getDomainMap());
    vec_type z(A.getRangeMap());
    vec_type resid(A.getRangeMap());

    z.randomize();

    // lambda: Current approximation of the eigenvalue of maximum magnitude.
    // normz: 2-norm of the current iteration vector z.
    // residual: 2-norm of the current residual vector 'resid'.
    //
    // Teuchos::ScalarTraits defines what zero and one means for any
    // type.  Most number types T know how to turn a 0 or a 1 (int)
    // into a T.  I have encountered some number types in C++ that do
    // not.  These tend to be extended-precision types that define
    // number operators and know how to convert from a float or
    // double, but don't have conversion operators for int.  Thus,
    // using Teuchos::ScalarTraits makes this code maximally general.
    scalar_type lambda = STS::zero ();
    magnitude_type normz = STM::zero ();
    magnitude_type residual = STM::zero ();

    const scalar_type one = STS::one();
    const scalar_type zero = STS::zero();

    const int reportFreq = 10;

    for(int iter = 0; iter < niters; iter++)
    {
      normz = z.norm2();
      q.scale(one/normz, z); // q := z/normz
      A.apply(q, z); // z := A*q

      lambda = q.dot(z);

      if(iter % reportFreq == 0 || iter+1 == niters)
      {
        resid.update(one, z, -lambda, q, zero);
        residual = resid.norm2();
        out << "iteration: " << iter;
        out << "\nlambda: " << lambda;
        out <<"||A*q - lamda*a||_2: " << residual << std::endl;
      }

      if(residual < tolerance)
      {
        out << "Converged after "<< iter << " iterations" << std::endl;
        break;
      }else if (iter +1 == niters)
      {
        out << "Failed to converge after "<< iter << " iterations" << std::endl;
        break;
      }

    }
    return lambda;
  }
};




void PowerMethdExample( Tpetra::global_size_t global_N,const RCP<const Teuchos::Comm<int>> &comm,std::ostream &out)
{

  // calculate dominate eigen mode via power method
  out <<"Setting " << global_N << "x" << global_N << " elements" << std::endl;
  // construct a map that puts the same number
  // of equations on each processor
  const global_ordinal_type idxBase = 0;
  RCP<const map_type> map = rcp(new map_type(global_N,idxBase,comm));

  const size_t myEls = map->getNodeNumElements();
  out << "Local els: " << myEls << std::endl;
  
  // lets see our local set of global indices
  ArrayView<const global_ordinal_type> myGlobalEls = map->getNodeElementList();

  comm->barrier();

//  printf("proc %d global id:{", comm->getRank());
//  for(auto &&gid : myGlobalEls)
//    printf("%d ", gid);
//  printf("}\n");

  // lets make a matrix!
  comm->barrier();
  out << "Create a sparse matrix!" << std::endl;

  // matrix has row distribution given by map
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type(map,0));

  // fill the sparse matrix one row at a time
  for(local_ordinal_type i = 0; i < (local_ordinal_type)myEls; i++)
  {
    // get global id for this local id
    const global_ordinal_type gid = map->getGlobalElement(i);
    if(gid == 0)
    {
      //A(0,0:10 = [2, -1]
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid, gid+1),
                            Teuchos::tuple<scalar_type>(2,-1));
    }else if(gid == global_N-1)
    {
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid-1,gid),
                            Teuchos::tuple<scalar_type>(-1,2));
    }else{
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid-1,gid, gid+1),
                            Teuchos::tuple<scalar_type>(-1,2,-1));
    }
  }

  // tell it that we are done adding etries to it
  A->fillComplete();

  // what are we dealing with
//  for(int i = 0; i < A->getNodeNumRows(); i++)
//  {
//    std::vector<local_ordinal_type> tmp;
//    ArrayView<const int> idxs;
//    ArrayView<const scalar_type> vals;
//    A->getLocalRowView(i, idxs, vals);
//
//    size_t c = 0;
//    for(auto idx : idxs)
//    {
//      printf("[%d, %d, %f]\n", i, idx, vals[c++]);
//    }
//  }

  // get estimate of dominate eigen value
  scalar_type lamb = PowerMethod<crs_matrix_type>::run(*A,500, 1e-2, out);
  out << "Lambda: " << lamb << std::endl;
}

void matrixExample( Tpetra::global_size_t global_N,const RCP<const Teuchos::Comm<int>> &comm, std::ostream &out)
{
  out <<"Setting " << global_N << "x" << global_N << " elements" << std::endl;
  // construct a map that puts the same number
  // of equations on each processor
  const global_ordinal_type idxBase = 0;
  RCP<const map_type> map = rcp(new map_type(global_N,idxBase,comm));

  const size_t myEls = map->getNodeNumElements();
  out << "Local els: " << myEls << std::endl;

  // lets see our local set of global indices
  ArrayView<const global_ordinal_type> myGlobalEls = map->getNodeElementList();

  comm->barrier();

  //  printf("proc %d global id:{", comm->getRank());
  //  for(auto &&gid : myGlobalEls)
  //    printf("%d ", gid);
  //  printf("}\n");

  // lets make a matrix!
  comm->barrier();
  out << "Create a sparse matrix!" << std::endl;
  
  // matrix has row distribution given by map
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type(map,0));

  // fill the sparse matrix one row at a time
  for(local_ordinal_type i = 0; i < (local_ordinal_type)myEls; i++)
  {
    // get global id for this local id
    const global_ordinal_type gid = map->getGlobalElement(i);
    if(gid == 0)
    {
      //A(0,0:10 = [2, -1]
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid, gid+1),
                            Teuchos::tuple<scalar_type>(1,1));
    }else if(gid == global_N-1)
    {
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid-1,gid),
                            Teuchos::tuple<scalar_type>(1,1));
    }else{
      A->insertGlobalValues(gid, Teuchos::tuple<global_ordinal_type>(gid-1,gid, gid+1),
                            Teuchos::tuple<scalar_type>(1,1,1));
    }
  }

  // tell it that we are done adding etries to it
  A->fillComplete();
  
  out << "A*A'" << std::endl;
  RCP<crs_matrix_type> B = rcp(new crs_matrix_type(map,0));
  Tpetra::MatrixMatrix::Multiply(*A, false, *A, true, *B);

   //what are we dealing with
    for(int i = 0; i < B->getNodeNumRows(); i++)
    {
      std::vector<local_ordinal_type> tmp;
      ArrayView<const int> idxs;
      ArrayView<const scalar_type> vals;
      B->getLocalRowView(i, idxs, vals);

      size_t c = 0;
      for(auto idx : idxs)
      {
        printf("B:[%d, %d, %d]\n", i, idx, vals[c++]);
      }
    }
  
  out << "\n B = A* A'" << std::endl;
  
  RCP<crs_graph_type> gB = rcp(new crs_graph_type(map,global_N));
  for(global_ordinal_type gid : myGlobalEls)
  {
    size_t numEntriesInRow = B->getNumEntriesInGlobalRow (gid);
    Array<scalar_type>         rowvals (numEntriesInRow);
    Array<global_ordinal_type> rowinds (numEntriesInRow);
    B->getGlobalRowCopy (gid, rowinds (), rowvals (), numEntriesInRow);
    printf("row %d = {",gid);
    for (size_t i = 0; i < numEntriesInRow; i++) {
      printf("%d ", rowvals[i]);
      if (rowvals[i] == 2)
        gB->insertGlobalIndices(gid, Teuchos::tuple<local_ordinal_type>(static_cast<local_ordinal_type>(rowinds[i])));
    }
    printf("}\n");
  }

  gB->fillComplete();
  out << "\nModifed B" << std::endl;

  out << "\nGet graph..." << std::endl;
  ArrayView<const global_ordinal_type> ggid = gB->getRowMap()->getNodeElementList();

  for(global_ordinal_type gid : ggid)
  {
    size_t numEntriesInRow = gB->getNumEntriesInGlobalRow(gid);
    Array<scalar_type>         rowvals (numEntriesInRow);
    Array<global_ordinal_type> rowinds (numEntriesInRow);
    gB->getGlobalRowCopy(gid,  rowinds, numEntriesInRow);

    printf("Neighbors %d = {",gid);
    for (size_t i = 0; i < numEntriesInRow; i++) {
      printf("%d ", rowinds[i]);
    }
    printf("}\n");
  }
}

int main(int narg, char ** arg)
{
  // Initialize Tpetra and get the default communicator.  Tpetra works
  // even if it was not built with MPI.  In that case, "comm" will be
  // some object that behaves like MPI_COMM_SELF.
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  Teuchos::oblackholestream blackHole;
  
  // get my rank and set a reference to the output stream
  // all processors not 0 speak into the black hole!
  const int myRank = comm->getRank ();
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // print current version of tpetra
  out << Tpetra::version() << std::endl;
  
  // do some shit
  Tpetra::global_size_t global_elements = 5;
  if(narg > 1) global_elements = atol(arg[1]);
  matrixExample(global_elements, comm, out);
  
  out << "FINISHED ROUTINE\n" << std::endl;
  return 0;
}
