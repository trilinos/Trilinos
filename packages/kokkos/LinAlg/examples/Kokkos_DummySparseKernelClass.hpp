#ifndef KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP
#define KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP

#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_BLAS_types.hpp>

/** \file Kokkos_DummySparseKernelClass.hpp
    \brief A file containing a stub for a new sparse kernel provider, as outlined in the \ref kokkos_crs_ops "Kokkos CRS API".
 */

namespace KokkosExamples {

  /** \class DummySparseKernel 
      \ingroup kokkos_crs_ops
      \brief A dummy-class illustrating the components necessary for a %Kokkos sparse operations provider.
   */
  template <class Node> 
  class DummySparseKernel {
  public:
    //@{ 
    //! @name Typedefs and structs

    //!
    typedef void  ScalarType;
    //!
    typedef void OrdinalType;
    //!
    typedef Node    NodeType;
    //! 
    typedef DummySparseKernel<Node> ThisType;

    /** \brief Rebind struct, for specifying type information for a different scalar.
        
        For a scalar type \c T, specify the class type that should be used for
        sparse kernels over sparse matrices of type \c T.
      */
    template <class T>
    struct rebind {
      typedef DummySparseKernel<Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! \brief Constructor accepting and retaining a node object.
    DummySparseKernel(const Teuchos::RCP<Node> & node) : node_(node) {}

    //! \brief Destructor.
    ~DummySparseKernel() {}

    //@} 
    //! @name Accessor routines.
    //@{ 

    /** \brief Node accessor.
        
        This returns the node type of type this::NodeType that was passed at object construction.
      */
    Teuchos::RCP<Node> getNode() const {return node_;}

    //@}
    //! @name Initialization of structure
    //@{

    /** \brief Initialize the structure of the sparse matrix.
        
        This is the mechanism by which the user specifies the structure for the sparse matrix. It always comes via 
        a Kokkos::CrsGraph<O,N,SO>, where
        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        The graph structure must be provided via initializeStructure() before the matrix values are specified to 
        initializeValues(). After calling initializeStructure(), the clear() method must be called before calling
        initializeStructure() again.

        In general, both initializeStructure() and initializeValues() must be called before calling multiply() or solve().
      */
    template <class Ordinal>
    void initializeStructure(const Kokkos::CrsGraph<Ordinal,Node,DummySparseKernel<Node> > & /* graph */) {}

    /** \brief Initialize the values of the sparse matrix.

        This is the mechanism by which the user specifies the values for the sparse matrix. It always comes via 
        a Kokkos::CrsMatrix<S,O,N,SO>, where
        - \c S is the scalar type this::ScalarType
        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        The graph structure must be provided via initializeStructure() before the matrix values are specified to 
        initializeValues(). The matrix values can be provided repeatadly via initializeValues() without the need
        to call clear().

        In general, both initializeStructure() and initializeValues() must be called before calling multiply() or solve().
      */
    template <class Scalar, class Ordinal>
    void initializeValues(const Kokkos::CrsMatrix<Scalar,Ordinal,Node,DummySparseKernel<Node> > & /* matrix */) {}

    /** \brief Clear the graph and matrix data.
        
        clear() must be called between successive calls to initializeStructure(). 

        In general, after calling clear(), no significant data (including persisting references) will be preserved, save for the pointer to the node. 

        In general, after calling clear(), multiply() and solve() cannot be called.
      */
    void clear() {}

    //@}
    //! @name Computational methods
    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const Kokkos::MultiVector<DomainScalar,Node> &X, Kokkos::MultiVector<RangeScalar,Node> &Y) const {}

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const Kokkos::MultiVector<DomainScalar,Node> &X, RangeScalar beta, Kokkos::MultiVector<RangeScalar,Node> &Y) const {}

    //! Solves the matrix for a given set of right-hand-sides.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const Kokkos::MultiVector<DomainScalar,Node> &Y, Kokkos::MultiVector<RangeScalar,Node> &X) const {}

    //@}
  protected:
    //! Copy constructor (protected and unimplemented)
    DummySparseKernel(const DummySparseKernel& source);

    Teuchos::RCP<Node> node_;
  };

  /** \example DummySparseKernelDriver.cpp 
    * This is an example that demonstrates the basic use case for a sparse kernel provider. It also verifies that the stub builds correctly.
    */

} // end of namespace KokkosExamples

namespace Kokkos {

  /** \brief A partial specialization of CrsGraph for use with KokkosExamples::DummySparseKernel
      \ingroup kokkos_crs_ops

      This specialization inherits from CrsGraphHostCompute. The consequence of this is that it is only appropriate for use on host-based nodes; and that 
      it doesn't specialize any capability based on the KokkosExamples::DummySparseKernel provider (which is fine, because this is for demonstration purposes only.)
    */
  template <class Ordinal, 
            class Node>
  class CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > : public CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > {
  public:
    CrsGraph(size_t numRows, const Teuchos::RCP<Node> &node) : CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &graph); // not implemented
  };

  /** \brief A partial specialization of CrsMatrix for use with KokkosExamples::DummySparseKernel
      \ingroup kokkos_crs_ops

      This specialization inherits from CrsMatrixHostCompute. The consequence of this is that it is only appropriate for use on host-based nodes; and that 
      it doesn't specialize any capability based on the KokkosExamples::DummySparseKernel provider (which is fine, because this is for demonstration purposes only.)
    */
  template <class Scalar,
            class Ordinal, 
            class Node>
  class CrsMatrix<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > : public CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > {
  public:
    CrsMatrix() : CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >() {}
    CrsMatrix(CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &graph) {
      const CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &hgraph = dynamic_cast<const CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >::setStaticGraph(hgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &graph) {
      CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &hgraph = dynamic_cast<CrsGraphHostCompute<Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> >::setOwnedGraph(hgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,Node,KokkosExamples::DummySparseKernel<Node> > &mat); // not implemented
  };

}

#endif
