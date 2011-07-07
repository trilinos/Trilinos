#ifndef MULTIPRECCG_HPP_
#define MULTIPRECCG_HPP_

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Kokkos_DefaultNode.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <iostream>
#include <functional>

#ifdef HAVE_TPETRA_QD
# include <qd/qd_real.h>
#endif

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

namespace TpetraExamples {

  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::Time;
  using Teuchos::null;
  using std::binary_function;
  using std::pair;
  using std::make_pair;
  using std::plus;
  using std::multiplies;

  struct trivial_fpu_fix {
    void fix() {}
    void unfix() {}
  };
  struct nontrivial_fpu_fix {
    unsigned int old_cw;
    void fix()   {fpu_fix_start(&old_cw);}
    void unfix() {fpu_fix_end(&old_cw);}
  };
  // implementations
  template <class T> struct fpu_fix : trivial_fpu_fix {};
#ifdef HAVE_TPETRA_QD
  template <> struct fpu_fix<qd_real> : nontrivial_fpu_fix {};
  template <> struct fpu_fix<dd_real> : nontrivial_fpu_fix {};
#endif

  //! Helper function to call CrsMatrix<T1>::convert<T2> if T1 != T2, and return the original matrix otherwise.
  template <class Tout, class Tin, class LO, class GO, class Node> 
  struct convertHelp {
    static RCP<const Tpetra::CrsMatrix<Tout,LO,GO,Node>> doit(const RCP<const Tpetra::CrsMatrix<Tin,LO,GO,Node>> &A)
    {
      return A->template convert<Tout>();
    }
  };

  //! Helper function to call CrsMatrix<T1>::convert<T2> if T1 != T2, and return the original matrix otherwise.
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

  //! \brief Non-member constructor to create a \c pair_op.
  template <class T1, class T2, class Op>
  pair_op<T1,T2,Op> make_pair_op(Op op) { return pair_op<T1,T2,Op>(op); }

  //! \brief Initialization for a database stack for the recursiveFPCG() algorithm, to be used with Tpetra::Ext::initStackDB().
  template <class S, class LO, class GO, class Node>
  class RFPCGInit 
  {
    private:
    typedef Tpetra::Map<LO,GO,Node>               Map;
    typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
    RCP<const CrsMatrix> A;

    public:

    RFPCGInit(RCP<Tpetra::CrsMatrix<S,LO,GO,Node>> Atop) : A(Atop) {}

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
      //
      db->set<RCP<const Op>>("A", AT                    );
      db->set("numIters", params.get<int>("numIters",A->getGlobalNumRows()) );
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
      ff.unfix();
      return db;
    }
  };

  /******************************
  *   Somewhat flexible CG
  *   Golub and Ye, 1999
  *
  *   r = b
  *   z = M*r
  *   p = z
  *   do
  *     alpha = r'*z / p'*A*p
  *     x = x + alpha*p
  *     r = r - alpha*A*p
  *     if outermost, check r for convergence
  *     z = M*r
  *     beta = z'*(r_new - r_old) / z'*r
  *     p = z + beta*p
  *   enddo
  ******************************/

  //! \brief Recursive, self-preconditioning flexible CG.
  template <class TS, class LO, class GO, class Node>      
  void recursiveFPCG(const RCP<Teuchos::FancyOStream> &out, ParameterList &db)
  {
    using Teuchos::as;
    using Tpetra::RTI::ZeroOp;
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
    RCP<const VectorT1> diag;
    if (TS::bottom) {
      diag = db.get<RCP<VectorT1>>("diag");
    }
    const T tolerance = db.get<double>("tolerance", 0.0);
    const int verbose = db.get<int>("verbose",0);
    static RCP<Time> timer;
    if (timer == null) {
      timer = Teuchos::TimeMonitor::getNewTimer(
                      "recurisveFPCG<"+Teuchos::TypeNameTraits<T>::name()+">"
              );
    }

    fpu_fix<T> ff;
    ff.fix();

    Teuchos::OSTab tab(out);

    if (verbose) {
      *out << "Beginning recursiveFPCG<" << Teuchos::TypeNameTraits<T>::name() << ">" << std::endl;
    }

    timer->start(); 
    const T r2 = TPETRA_BINARY_PRETRANSFORM_REDUCE(
                    r, x,                                         // fused: 
                    x,                                            //      : r = x  
                    r*r, ZeroOp<T>, plus<T>() );                  //      : sum r'*r
    // b comes in, x goes out. now we're done with b, so zero the solution.
    TPETRA_UNARY_TRANSFORM( x,  ST::zero() );                     // x = 0
    if (TS::bottom) {
      TPETRA_TERTIARY_TRANSFORM( z, diag, r,    r/diag );         // z = D\r
    }
    else {
      ff.unfix();
      timer->stop(); 
      ParameterList &db_T2 = db.sublist("child");
      auto bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
      TPETRA_BINARY_TRANSFORM( bx_T2, r, as<T2>(r) );             // b_T2 = (T2)r
      recursiveFPCG<typename TS::next,LO,GO,Node>(out,db_T2);     // x_T2 = A_T2 \ b_T2 
      TPETRA_BINARY_TRANSFORM( z, bx_T2, as<T>(bx_T2) );          // z    = (T)bx_T2
      timer->start(); 
      ff.fix();
    }
    T zr = TPETRA_TERTIARY_PRETRANSFORM_REDUCE( 
                    p, z, r,                                      // fused: 
                    z,                                            //      : p = z
                    z*r, ZeroOp<T>, plus<T>() );                  //      : z*r
    ///////////////////////////////////
    int k;
    for (k=0; k<numIters; ++k) 
    {
      A->apply(*p,*Ap);                                           // Ap = A*p
      T pAp = TPETRA_REDUCE2( p, Ap,     
                              p*Ap, ZeroOp<T>, plus<T>() );       // p'*Ap
      const T alpha = zr / pAp;
      TPETRA_BINARY_TRANSFORM( x,    p,  x + alpha*p  );          // x = x + alpha*p
      TPETRA_BINARY_TRANSFORM( rold, r,  r            );          // rold = r
      T rr = TPETRA_BINARY_PRETRANSFORM_REDUCE(
                               r, Ap,                             // fused:
                               r - alpha*Ap,                      //      : r - alpha*Ap
                               r*r, ZeroOp<T>, plus<T>() );       //      : sum r'*r
      if (verbose > 1) *out << "|res|/|res_0|: " << ST::squareroot(rr/r2) 
                            << std::endl;
      if (rr/r2 < tolerance*tolerance) {
        if (verbose) {
          *out << "Convergence detected!" << std::endl;
        }
        break;
      }
      if (TS::bottom) {
        TPETRA_TERTIARY_TRANSFORM( z, diag, r,    r/diag );        // z = D\r
      }
      else {
        ff.unfix();
        timer->stop(); 
        ParameterList &db_T2 = db.sublist("child");
        auto bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
        TPETRA_BINARY_TRANSFORM( bx_T2, r,    as<T2>(r) );         // b_T2 = (T2)r
        recursiveFPCG<typename TS::next,LO,GO,Node>(out,db_T2);    // x_T2 = A_T2 \ b_T2
        TPETRA_BINARY_TRANSFORM( z, bx_T2,    as<T>(bx_T2) );      // z    = (T)bx_T2
        timer->start(); 
        ff.fix();
      }
      const T zoro = zr;                                                         
      typedef ZeroOp<pair<T,T>> ZeroPTT;
      auto plusTT = make_pair_op<T,T>(plus<T>());
      pair<T,T> both = TPETRA_REDUCE3( z, r, rold,                 // fused: z'*r and z'*r_old
                                       make_pair(z*r, z*rold), ZeroPTT, plusTT );
      zr = both.first; // this is used on the next iteration as well
      const T znro = both.second;
      const T beta = (zr - znro) / zoro;
      TPETRA_BINARY_TRANSFORM( p, z,   z + beta*p );               // p = z + beta*p
    }
    timer->stop(); 
    ff.unfix();
    if (verbose) {
      *out << "Leaving recursiveFPCG<" << Teuchos::TypeNameTraits<T>::name() 
           << "> after " << k << " iterations." << std::endl;
    }
  }

  //! \brief Recursive, self-preconditioning flexible CG.
  template <class TS, class LO, class GO, class Node>      
  void recursiveFPCGUnfused(const RCP<Teuchos::FancyOStream> &out, ParameterList &db)
  {
    using Tpetra::RTI::unary_transform;
    using Tpetra::RTI::binary_transform;
    using Tpetra::RTI::tertiary_transform;
    using Tpetra::RTI::reduce;
    using Tpetra::RTI::reductionGlob;
    using Tpetra::RTI::ZeroOp;
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
    static RCP<Time> timer;
    if (timer == null) {
      timer = Teuchos::TimeMonitor::getNewTimer(
                      "recursiveFPCGUnfused<"+Teuchos::TypeNameTraits<T>::name()+">"
              );
    }
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
      *out << "Beginning recursiveFPCGUnfused<" << Teuchos::TypeNameTraits<T>::name() << ">" << std::endl;
    }

    timer->start(); 
    binary_transform( *r, *x,            [](T, T bi)             {return bi;});  // r = b     (b is stored in x)
    const T r2 = reduce( *r, *r,      reductionGlob<ZeroOp<T>>(multiplies<T>(),  // r'*r
                                                                   plus<T>())); 
    unary_transform(  *x,                [](T)           {return ST::zero();});  // set x = 0 (now that we don't need b)
    if (TS::bottom) {
      tertiary_transform( *z, *diag, *r, [](T, T di, T ri)    {return ri/di;});  // z = D\r
    }
    else {
      ff.unfix();
      timer->stop(); 
      ParameterList &db_T2 = db.sublist("child");
      auto bx_T2 = db_T2.get< RCP<VectorT2>>("bx");
      binary_transform( *bx_T2, *r, [](T2, T ri)         {return as<T2>(ri);});  // b_T2 = (T2)r       
      recursiveFPCGUnfused<typename TS::next,LO,GO,Node>(out,db_T2);             // x_T2 = A_T2 \ b_T2 
      binary_transform( *z, *bx_T2, [](T, T2 xi)          {return as<T>(xi);});  // z    = (T)x_T2     
      timer->start(); 
      ff.fix();
    }
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
        recursiveFPCGUnfused<typename TS::next,LO,GO,Node>(out,db_T2);           // x_T2 = A_T2 \ b_T2
        binary_transform( *z, *bx_T2, [](T, T2 xi)        {return as<T>(xi);});  // z    = (T)x_T2
        timer->start(); 
        ff.fix();
      }
      const T zoro = zr;                                                         
      zr = reduce( *z, *r,          reductionGlob<ZeroOp<T>>(multiplies<T>(),    // z'*r
                                                             plus<T>()));        // this is loop-carried
      const T znro = reduce( *z, *rold, reductionGlob<ZeroOp<T>>(multiplies<T>(),// z'*r_old
                                                                 plus<T>())); 
      const T beta = (zr - znro) / zoro;
      binary_transform( *p, *z, [beta](T pi, T zi)     {return zi + beta*pi;});  // p = z + beta*p
    }
    timer->stop(); 
    ff.unfix();
    if (verbose) {
      *out << "Leaving recursiveFPCGUnfused<" << Teuchos::TypeNameTraits<T>::name() << "> after " << k << " iterations." << std::endl;
    }
  }

} // end of namespace TpetraExamples

#endif // MULTIPRECCG_HPP_
