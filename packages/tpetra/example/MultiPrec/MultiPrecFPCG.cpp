#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Kokkos_DefaultNode.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_TypeStack.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <iostream>
#include <iomanip>
#include <functional>

#ifdef HAVE_TPETRA_QD
#include <qd/qd_real.h>
#endif

/** \file MultiPrecFPCG.cpp
    \brief An example of a multi-precision algorithm, using a flexible preconditioned CG with recursive precision.
 */

namespace Tpetra {
  namespace RTI {
    // specialization for pair
    template <class T1, class T2>
    class ZeroOp<std::pair<T1,T2>> {
      public:
      static inline std::pair<T1,T2> identity() {
        return std::make_pair( Teuchos::ScalarTraits<T1>::zero(), 
                               Teuchos::ScalarTraits<T2>::zero() );
      }
    };
  }
}

namespace RFPCG {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Tpetra::RTI::unary_transform;
    using Tpetra::RTI::binary_transform;
    using Tpetra::RTI::reduce;
    using Tpetra::RTI::reductionGlob;
    using Tpetra::RTI::ZeroOp;
    using std::binary_function;
    using std::pair;
    using std::make_pair;
    using std::plus;
    using std::multiplies;
  /** \class pair_op 
    \brief pair_op is a reduction function object that takes two arguments of a specified precision and multiplies them using a different precision.
   */
  template <class T1, class T2, class Op>
  class pair_op : public binary_function<pair<T1,T2>,pair<T1,T2>,pair<T1,T2>> {
  private:
    Op op_;
  public:
    pair_op(Op op) : op_(op) {}
    inline pair<T1,T2> operator()(const pair<T1,T2>& a, const pair<T1,T2>& b) const {
      return make_pair(op_(a.first,b.first),op_(a.second,b.second));
    }
  };

  template <class T1, class T2, class Op>
  pair_op<T1,T2,Op> make_pair_op(Op op) { return pair_op<T1,T2,Op>(op); }

  enum IterationLevel {
    OuterMostLevel,
    NotOuterMostLevel
  };

  template <class S, class LO, class GO, class Node>
  class CGInit 
  {
    private:
    typedef Tpetra::Map<LO,GO,Node>               Map;
    typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
    RCP<const CrsMatrix> A;
  
    public:
  
    CGInit(RCP<Tpetra::CrsMatrix<S,LO,GO,Node>> Atop) : A(Atop) {}
  
    template <class T>
    RCP<ParameterList> initDB(ParameterList &params) 
    {
      typedef Tpetra::Vector<T,LO,GO,Node>    Vector;
      typedef Tpetra::Operator<T,LO,GO,Node>      Op;
      RCP<const Map> map = A->getDomainMap();
      //
      RCP<ParameterList> db = Teuchos::parameterList();
      db->set< int         >("numIters", params.get<int>("numIters",5) );
      db->set<RCP<Vector>>("bx",   Tpetra::createVector<T>(map) );
      db->set<RCP<Vector>>("r",    Tpetra::createVector<T>(map) );
      db->set<RCP<Vector>>("z",    Tpetra::createVector<T>(map) );
      db->set<RCP<Vector>>("p",    Tpetra::createVector<T>(map) );
      db->set<RCP<Vector>>("Ap",   Tpetra::createVector<T>(map) );
      db->set<RCP<Vector>>("rold", Tpetra::createVector<T>(map) );
      db->set<RCP<const Op>>("A",  A->template convert<T>()     );
      return db;
    }
  };

  template <class TS, class LO, class GO, class Node>      
  void recursiveFPCG(IterationLevel ilevel, ParameterList &db)
  {
    typedef typename TS::type       T;
    typedef typename TS::next::type T2;
    typedef Tpetra::Vector<T,LO,GO,Node> VectorT1;
    typedef Tpetra::Vector<T2,LO,GO,Node> VectorT2;
    typedef Tpetra::Operator<T,LO,GO,Node>   OpT1;
    typedef Teuchos::ScalarTraits<T>            ST;
    // get objects from level database
    const int numIters = db.get<int>("numIters");
    RCP<VectorT1> r    = db.get<RCP<VectorT1>>("r");
    RCP<VectorT1> z    = db.get<RCP<VectorT1>>("z");
    RCP<VectorT1> p    = db.get<RCP<VectorT1>>("p");
    RCP<VectorT1> Ap   = db.get<RCP<VectorT1>>("Ap");
    RCP<VectorT1> rold = db.get<RCP<VectorT1>>("rold");
    RCP<VectorT1> x    = db.get<RCP<VectorT1>>("bx");
    RCP<const OpT1> A  = db.get<RCP<const OpT1>>("A");
    const T epsilon = db.get<T>("tol", ST::zero());
  
    /* 
        r = b
        z = M*r
        p = z
        do
          alpha = r'*z / p'*A*p
          x = x + alpha*p
          r = r - alpha*A*p
          if outermost, check r for convergence
          z = M*r
          beta = z'*(r_new - r_old) / z'*r
          p = z + beta*p
        enddo
     */
    binary_transform( *r, *x, [](T, T bi)                        {return bi;});  // r = b
    unary_transform( *x,      [](T)                      {return ST::zero();});  // set x = 0 (now that we don't need b)
    if (TS::bottom) {
      binary_transform( *z, *r, [](T, T ri)                      {return ri;});  // z = r
    }
    else {
      ParameterList &db_T2 = db.sublist("child");
      RCP<VectorT2> bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
      binary_transform( *bx_T2, *r, [](T2, T ri)                 {return ri;});  // b_T2 = r
      recursiveFPCG<typename TS::next,LO,GO,Node>(NotOuterMostLevel,db_T2);      // x_T2 = A_T2 \ b_T2
      binary_transform( *z, *bx_T2, [](T, T2 xi)                 {return xi;});  // z    = x_T2
    }
    binary_transform( *p, *z, [](T, T zi)                        {return zi;});  // p = z
    T zr = reduce( *z, *r,   reductionGlob<ZeroOp<T>>(multiplies<T>(),           // z'*r
                                                            plus<T>())); 
    ///////////////////////////////////
    for (int k=0; k<numIters; ++k) 
    {
      A->apply(*p,*Ap);                                                          // Ap = A*p
      T pAp = reduce( *p, *Ap, reductionGlob<ZeroOp<T>>(multiplies<T>(),         // p'*Ap
                                                              plus<T>())); 
      const T alpha = zr / pAp;
      binary_transform( *x, *p, [alpha](T xi, T pi)   {return xi + alpha*pi;});  // x = x + alpha*p
      binary_transform( *rold, *r, [](T, T ri)                   {return ri;});  // rold = r
      binary_transform( *r, *Ap,[alpha](T ri, T Api) {return ri - alpha*Api;});  // r = r - alpha*Ap
      if (ilevel == OuterMostLevel) {
        T rr = reduce( *r, *r, reductionGlob<ZeroOp<T>>(multiplies<T>(),         // r'*r
                                                              plus<T>())); 
        if (rr < epsilon*epsilon) break;
      }
      if (TS::bottom) {
        binary_transform( *z, *r, [](T, T ri)                    {return ri;});  // z = r
      }
      else {
        ParameterList &db_T2 = db.sublist("child");
        RCP<VectorT2> bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
        binary_transform( *bx_T2, *r, [](T2, T ri)               {return ri;});  // b_T2 = r
        recursiveFPCG<typename TS::next,LO,GO,Node>(NotOuterMostLevel,db_T2);    // x_T2 = A_T2 \ b_T2
        binary_transform( *z, *bx_T2, [](T, T2 xi)               {return xi;});  // z    = x_T2
      }
      const T zoro = zr;                                                         // simultaneous
      pair<T,T> both = reduce(*z, *r, *rold,                                     // z'*r     and
                           reductionGlob<ZeroOp<pair<T,T>>>(                     // z'*r_old
                                         [](T zi, T ri, T roldi) { return make_pair(zi*ri, zi*roldi); },
                                         make_pair_op<T,T>(plus<T>()))
                    );
      zr = both.first; // this is used on the next iteration as well
      const T znro = both.second;
      const T beta = (zr - znro) / zoro;
      binary_transform( *p, *z, [beta](T pi, T zi)     {return zi + beta*pi;});  // p = z + beta*p
    }
  }
}

int main(int argc, char *argv[])
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  using Teuchos::TimeMonitor;
  using Teuchos::RCP;
  using Teuchos::ParameterList;

// #ifdef HAVE_TPETRA_QD
//   TPETRAEXT_TYPESTACK4(MPStack, qd_real, dd_real, double, float )
// #else
  TPETRAEXT_TYPESTACK2(MPStack, double, float )
// #endif

  // 
  // Get the default communicator and node
  //
  auto &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  auto comm = platform.getComm();
  auto node = platform.getNode();
  const int myImageID = comm->getRank();

  // Static types
  typedef typename MPStack::type   S;
  typedef int                     LO;
  typedef int                     GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  std::string matfile;
  std::string xmlfile;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("matrix-file",&matfile,"Filename for matrix");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl;
    std::cout << "Comm info: " << *comm;
    std::cout << "Node type: " << Teuchos::typeName(*node) << std::endl;
    std::cout << "Problem scalar: " << Teuchos::TypeNameTraits<S>::name() << std::endl;
    std::cout << std::endl;
  }

  // read the matrix
  RCP<CrsMatrix> A;
  Tpetra::Utils::readHBMatrix(matfile,comm,node,A);

  // get the solver parameters
  Teuchos::ParameterList stackPL;
  std::string xmlString(
    " <ParameterList>                                          \n"
    " <Parameter name='numIters' value='5' type='int'/>        \n"
    "   <ParameterList name='child'>                           \n"
    "   <Parameter name='numIters' value='5' type='int'/>      \n"
    "     <ParameterList name='child'>                         \n"
    "     <Parameter name='numIters' value='5' type='int'/>    \n"
    "       <ParameterList name='child'>                       \n"
    "       <Parameter name='numIters' value='5' type='int'/>  \n"
    "       </ParameterList>                                   \n"
    "     </ParameterList>                                     \n"
    "   </ParameterList>                                       \n"
    " </ParameterList>                                         \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,&stackPL);
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,&stackPL);

  // init the solver stack
  RFPCG::CGInit<S,LO,GO,Node> init(A);
  RCP<ParameterList> db = Tpetra::Ext::initStackDB<MPStack>(stackPL,init);

  // call the solve
  RFPCG::recursiveFPCG<MPStack,LO,GO,Node>(RFPCG::OuterMostLevel,*db);

  //
  // Print timings
  //
  if (verbose) {
    TimeMonitor::summarize( std::cout );
  }

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}

/** \example MultiPrecFPCG.cpp
    Demonstrate using Tpetra::RTI and a multi-precision flexible preconditioned CG.
  */
