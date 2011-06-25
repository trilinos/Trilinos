#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>

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

struct trivial_fpu_fix {
  void fix() {}
  void unfix() {}
};

struct nontrivial_fpu_fix {
  unsigned int old_cw;
  void fix()   {fpu_fix_start(&old_cw);}
  void unfix() {fpu_fix_end(&old_cw);}
};

template <class T>
struct fpu_fix : trivial_fpu_fix {};

#ifdef HAVE_TPETRA_QD
template <>
struct fpu_fix<qd_real> : nontrivial_fpu_fix {};
template <>
struct fpu_fix<dd_real> : nontrivial_fpu_fix {};
#endif


namespace RFPCG {
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::Time;
  using Tpetra::RTI::unary_transform;
  using Tpetra::RTI::binary_transform;
  using Tpetra::RTI::tertiary_transform;
  using Tpetra::RTI::reduce;
  using Tpetra::RTI::reductionGlob;
  using Tpetra::RTI::ZeroOp;
  using std::binary_function;
  using std::pair;
  using std::make_pair;
  using std::plus;
  using std::multiplies;

  template <class Tout, class Tin, class LO, class GO, class Node> 
  struct convertHelp {
    static RCP<const Tpetra::CrsMatrix<Tout,LO,GO,Node>> doit(const RCP<const Tpetra::CrsMatrix<Tin,LO,GO,Node>> &A)
    {
      return A->template convert<Tout>();
    }
  };

  template <class T, class LO, class GO, class Node> 
  struct convertHelp<T,T,LO,GO,Node> {
    static RCP<const Tpetra::CrsMatrix<T,LO,GO,Node>> doit(const RCP<const Tpetra::CrsMatrix<T,LO,GO,Node>> &A)
    {
      return A;
    }
  };

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
      fpu_fix<T> ff;
      ff.fix();
      typedef Tpetra::Vector<T,LO,GO,Node>    Vector;
      typedef Tpetra::Operator<T,LO,GO,Node>      Op;
      typedef Tpetra::CrsMatrix<T,LO,GO,Node>    Mat;
      RCP<const Map> map = A->getDomainMap();
      RCP<ParameterList> db = Teuchos::parameterList();
      RCP<const Mat> AT = convertHelp<T,S,LO,GO,Node>::doit(A);
      auto timer = Teuchos::TimeMonitor::getNewTimer(
                      "recursiveFPCG<"+Teuchos::TypeNameTraits<T>::name()+">"
                   );
      //
      db->set<RCP<const Op>>("A", AT                    );
      db->set("numIters", params.get<int>("numIters",5) );
      db->set("tolerance",params.get<double>("tolerance",1e-7));
      db->set("verbose",  params.get<int>("verbose",0) );
      db->set("bx",       Tpetra::createVector<T>(map)  );
      db->set("r",        Tpetra::createVector<T>(map)  );
      db->set("z",        Tpetra::createVector<T>(map)  );
      db->set("p",        Tpetra::createVector<T>(map)  );
      db->set("Ap",       Tpetra::createVector<T>(map)  );
      db->set("rold",     Tpetra::createVector<T>(map)  );
      if (params.get<bool>("Extract Diagonal",false)) {
        RCP<Vector> diag = Tpetra::createVector<T>(map);
        AT->getLocalDiagCopy(*diag);
        db->set("diag", diag);
      }
      db->set("timer", timer);
      ff.unfix();
      return db;
    }
  };

  template <class TS, class LO, class GO, class Node>      
  void recursiveFPCG(const RCP<Teuchos::FancyOStream> &out, ParameterList &db)
  {
    using Teuchos::as;
    typedef typename TS::type       T;
    typedef typename TS::next::type T2;
    typedef Tpetra::Vector<T ,LO,GO,Node> VectorT1;
    typedef Tpetra::Vector<T2,LO,GO,Node> VectorT2;
    typedef Tpetra::Operator<T,LO,GO,Node>    OpT1;
    typedef Teuchos::ScalarTraits<T>            ST;
    // get objects from level database
    const int numIters = db.get<int>("numIters");
    auto x     = db.get<RCP<VectorT1>>("bx");
    auto r     = db.get<RCP<VectorT1>>("r");
    auto z     = db.get<RCP<VectorT1>>("z");
    auto p     = db.get<RCP<VectorT1>>("p");
    auto Ap    = db.get<RCP<VectorT1>>("Ap");
    auto rold  = db.get<RCP<VectorT1>>("rold");
    auto A     = db.get<RCP<const OpT1>>("A");
    auto timer = db.get<RCP<Time>>("timer");
    RCP<const VectorT1> diag;
    if (TS::bottom) {
      diag = db.get<RCP<VectorT1>>("diag");
    }
    const T tolerance = db.get<double>("tolerance", 0.0);
    const int verbose = db.get<int>("verbose",0);

    fpu_fix<T> ff;
    ff.fix();

    Teuchos::OSTab tab(out);

    if (verbose) {
      *out << "Beginning recursiveFPCG<" << Teuchos::TypeNameTraits<T>::name() << ">" << std::endl;
    }

    timer->start(); 
    /* 
        Somewhat flexible CG

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
    // FINISH: combine
    binary_transform( *r, *x, [](T, T bi)                        {return bi;});  // r = b     (b is stored in x)
    const T r2 = reduce( *r, *r,    reductionGlob<ZeroOp<T>>(multiplies<T>(),    // r'*r
                                                             plus<T>())); 
    unary_transform(  *x,     [](T)                      {return ST::zero();});  // set x = 0 (now that we don't need b)
    if (TS::bottom) {
      tertiary_transform( *z, *diag, *r, [](T, T di, T ri)    {return ri/di;});  // z = D\r
    }
    else {
      ff.unfix();
      timer->stop(); 
      ParameterList &db_T2 = db.sublist("child");
      auto bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
      binary_transform( *bx_T2, *r, [](T2, T ri)         {return as<T2>(ri);});  // b_T2 = (T2)r       
      recursiveFPCG<typename TS::next,LO,GO,Node>(out,db_T2);                    // x_T2 = A_T2 \ b_T2 
      binary_transform( *z, *bx_T2, [](T, T2 xi)          {return as<T>(xi);});  // z    = (T)x_T2     
      timer->start(); 
      ff.fix();
    }
    // FINISH: combine
    binary_transform( *p, *z, [](T, T zi)                        {return zi;});  // p = z
    T zr = reduce( *z, *r,          reductionGlob<ZeroOp<T>>(multiplies<T>(),    // z'*r
                                                                   plus<T>())); 
    ///////////////////////////////////
    int k;
    for (k=0; k<numIters; ++k) 
    {
      A->apply(*p,*Ap);                                                          // Ap = A*p
      T pAp = reduce( *p, *Ap,      reductionGlob<ZeroOp<T>>(multiplies<T>(),    // p'*Ap
                                                                   plus<T>())); 
      const T alpha = zr / pAp;
      binary_transform( *x, *p, [alpha](T xi, T pi)   {return xi + alpha*pi;});  // x = x + alpha*p
      binary_transform( *rold, *r, [](T, T ri)                   {return ri;});  // rold = r
      // FINISH: combine
      binary_transform( *r, *Ap,[alpha](T ri, T Api) {return ri - alpha*Api;});  // r = r - alpha*Ap
      T rr = reduce( *r, *r,      reductionGlob<ZeroOp<T>>(multiplies<T>(),      // r'*r
                                                                 plus<T>())); 
      if (verbose > 1) *out << "|res|/|res_0|: " << ST::squareroot(rr/r2) << std::endl;
      if (rr/r2 < tolerance*tolerance) {
        if (verbose) {
          *out << "Convergence detected!" << std::endl;
        }
        break;
      }
      if (TS::bottom) {
        tertiary_transform( *z, *diag, *r, [](T, T di, T ri)  {return ri/di;});  // z = D\r
      }
      else {
        ff.unfix();
        timer->stop(); 
        ParameterList &db_T2 = db.sublist("child");
        auto bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
        binary_transform( *bx_T2, *r, [](T2, T ri)       {return as<T2>(ri);});  // b_T2 = (T2)r
        recursiveFPCG<typename TS::next,LO,GO,Node>(out,db_T2);                  // x_T2 = A_T2 \ b_T2
        binary_transform( *z, *bx_T2, [](T, T2 xi)        {return as<T>(xi);});  // z    = (T)x_T2
        timer->start(); 
        ff.fix();
      }
      const T zoro = zr;                                                         // simultaneous
      pair<T,T> both = reduce(*z, *r, *rold,                                     // z'*r     and
                              reductionGlob<ZeroOp<pair<T,T>>>(                  // z'*r_old
                              [](T zi, T ri, T roldi) { return make_pair(zi*ri, zi*roldi); },
                              make_pair_op<T,T>(plus<T>())) );
      zr = both.first; // this is used on the next iteration as well
      const T znro = both.second;
      const T beta = (zr - znro) / zoro;
      binary_transform( *p, *z, [beta](T pi, T zi)     {return zi + beta*pi;});  // p = z + beta*p
    }
    timer->stop(); 
    ff.unfix();
    if (verbose) {
      *out << "Leaving recursiveFPCG<" << Teuchos::TypeNameTraits<T>::name() << "> after " << k << " iterations." << std::endl;
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
  using std::pair;
  using std::make_pair;
  using std::plus;
  using std::cout;
  using std::endl;
  using RFPCG::make_pair_op;
  using Tpetra::RTI::reductionGlob;
  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::binary_pre_transform_reduce;
  using Tpetra::RTI::binary_transform;

#ifdef HAVE_TPETRA_QD
  TPETRAEXT_TYPESTACK3(MPStack, qd_real, dd_real, double )
#else
  TPETRAEXT_TYPESTACK2(MPStack, double, float )
#endif

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
  typedef Tpetra::Vector<S,LO,GO,Node>       Vector;

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

  Teuchos::RCP<Teuchos::FancyOStream> out; 
  if (verbose) out = Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false));
  else         out = Teuchos::getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));

  // 
  // Say hello, print some communicator info
  //
  *out << "\n" << Tpetra::version() << endl
       << "Comm info: " << *comm
       << "Node type: " << Teuchos::typeName(*node) << endl
       << "Outermost scalar: " << Teuchos::TypeNameTraits<S>::name() << endl
       << endl;

  // read the matrix
  RCP<CrsMatrix> A;
  Tpetra::Utils::readHBMatrix(matfile,comm,node,A);

  // get the solver parameters
  Teuchos::ParameterList stackPL;
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                           \n"
    " <Parameter name='numIters' value='10' type='int'/>        \n"
    "   <ParameterList name='child'>                            \n"
    "   <Parameter name='numIters' value='20' type='int'/>      \n"
    "     <ParameterList name='child'>                          \n"
    "     <Parameter name='numIters' value='40' type='int'/>    \n"
    "       <ParameterList name='child'>                        \n"
    "       <Parameter name='numIters' value='80' type='int'/>  \n"
    "       </ParameterList>                                    \n"
    "     </ParameterList>                                      \n"
    "   </ParameterList>                                        \n"
    " </ParameterList>                                          \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,&stackPL);
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,&stackPL);

  // init the solver stack
  RFPCG::CGInit<S,LO,GO,Node> init(A);
  RCP<ParameterList> db = Tpetra::Ext::initStackDB<MPStack>(stackPL,init);
  
  // choose a solution, compute a right-hand-side
  auto x = Tpetra::createVector<S>(A->getRowMap()),
       b = Tpetra::createVector<S>(A->getRowMap());
  x->randomize();
  A->apply(*x,*b);
  {
    auto bx = db->get<RCP<Vector>>("bx");
    binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
  }

  // call the solve
  RFPCG::recursiveFPCG<MPStack,LO,GO,Node>(out,*db);

  //
  // Print timings
  //
  TimeMonitor::summarize( *out );

  //
  // Test the result
  //
  auto xhat = db->get<RCP<Vector>>("bx"),
       bhat = Tpetra::createVector<S>(A->getRowMap());
  A->apply(*xhat,*bhat);
  // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
  auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                    reductionGlob<ZeroOp<pair<S,S>>>( 
                                    [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                    [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                    make_pair_op<S,S>(plus<S>())) );
  const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
          bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
  // check that residual is as requested
  *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
  const double tolerance = db->get<double>("tolerance");
  *out << "tolerance: " << tolerance << endl;
  if (enrm / bnrm < tolerance) {
    *out << "End Result: TEST PASSED" << endl;
  }
  return 0;
}

/** \example MultiPrecFPCG.cpp
    Demonstrate using Tpetra::RTI and a multi-precision flexible preconditioned CG.
  */
