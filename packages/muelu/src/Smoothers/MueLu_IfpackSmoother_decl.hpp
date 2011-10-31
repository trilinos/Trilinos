#ifndef MUELU_IFPACKSMOOTHER_DECL_HPP
#define MUELU_IFPACKSMOOTHER_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK

#include <Ifpack.h>

#include <Epetra_CrsMatrix.h>

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class IfpackSmoother
    @brief Class that encapsulates Ifpack smoothers.
    
    This class creates an Ifpack preconditioner factory.  The factory creates a smoother based on the
    type and ParameterList passed into the constructor.  See the constructor for more information.
  */
  class IfpackSmoother : public SmootherPrototype<double,int,int>
  {

    typedef double Scalar;
    typedef int    LocalOrdinal;
    typedef int    GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

        The options passed into IfpackSmoother are those given in the Ifpack user's manual.

        @param type smoother type
        @param list options for the particular smoother (e.g., fill factor or damping parameter)

        Here is how to select some of the most common smoothers.

         - Gauss-Seidel
            - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
            - parameter list options
                - <tt>relaxation: type</tt> = <tt>Gauss-Seidel</tt>
                - <tt>relaxation: damping factor</tt>
         - symmetric Gauss-Seidel
            - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
            - parameter list options
                - <tt>relaxation: type</tt> = <tt>symmetric Gauss-Seidel</tt>
                - <tt>relaxation: damping factor</tt>
         - Chebyshev
            - <tt>type</tt> = <tt>Chebyshev</tt>
            - parameter list options
                - <tt>chebyshev: ratio eigenvalue</tt>
                - <tt>chebyshev: min eigenvalue</tt>
                - <tt>chebyshev: max eigenvalue</tt>
                - <tt>chebyshev: degree</tt>
                - <tt>chebyshev: zero starting solution</tt> (defaults to <tt>true</tt>)
         - ILU
            - <tt>type</tt> = <tt>ILU</tt>
            - parameter list options
                - <tt>fact: level-of-fill</tt>

        See also Ifpack_PointRelaxation, Ifpack_Chebyshev, Ifpack_ILU.
    */
    IfpackSmoother(std::string const & type, Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0, RCP<FactoryBase> AFact = Teuchos::null) //TODO: empty paramList valid for Ifpack??
      : type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
    ;

    //! Destructor
    virtual ~IfpackSmoother() ;

    //@}

    //! @name Set/Get methods

    //@{

    //! Set smoother parameters
    void SetParameters(Teuchos::ParameterList const & paramList) ;

    //! Get smoother parameters
    Teuchos::ParameterList const & GetParameters() ;

    //JG: I'm not sure if it's a good idea to provide Get/Set NIts (for code maintainability)
    
    //     /*! @brief Set the number of smoothing sweeps/degree.
    //
    //        If the smoother is relaxation, this sets the number of sweeps.
    //        If the smoother is Chebyshev, this sets the polynomial degree.
    //
    //        Note:  This can be called after the preconditioner is set up, i.e., after
    //        calling IfpackSmoother::Setup().
    //     */
    //     void SetNIts(LO const &nIts) ;
    //
    //     /*! @brief Get the number of smoothing sweeps.
    //
    //        If the smoother is relaxation, this returns the number of sweeps.
    //        If the smoother is Chebyshev, this returns the polynomial degree.
    //     */
    //     LO GetNIts() ;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

        This creates the underlying Ifpack smoother object, copies any parameter list options
        supplied to the constructor to the Ifpack object, and computes the preconditioner.
    */
    void Setup(Level &currentLevel) ;

    /*! @brief Apply the preconditioner.

        Solves the linear system <tt>AX=B</tt> using the constructed smoother.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero (optional) If false, some work can be avoided.  Whether this actually saves any work depends on the underlying Ifpack implementation.
    */
    void Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero=false) const ;

    //@}

    //! @name Utilities
    //@{

    RCP<SmootherPrototype> Copy() const ;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const ;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;; // class IfpackSmoother

  //! Non-member templated function GetIfpackSmoother() returns a new IfpackSmoother object when <Scalar, LocalOrdinal, GlobalOrdinal> == <double, int, int>. Otherwise, an exception is thrown.
  //! This function simplifies the usage of IfpackSmoother objects inside of templates as templates do not have to be specialized for <double, int, int> (see DirectSolver for an example).
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > GetIfpackSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LocalOrdinal const &overlap=0, RCP<FactoryBase> AFact = Teuchos::null) ;
  //
  template <>
  inline RCP<MueLu::SmootherPrototype<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > GetIfpackSmoother<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>(std::string const & type, Teuchos::ParameterList const & paramList, int const &overlap, RCP<FactoryBase> AFact) ;

} // namespace MueLu

#define MUELU_IFPACK_SMOOTHER_SHORT
#endif // ifdef HAVE_MUELU_IFPACK
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_IFPACKSMOOTHER_DECL_HPP
