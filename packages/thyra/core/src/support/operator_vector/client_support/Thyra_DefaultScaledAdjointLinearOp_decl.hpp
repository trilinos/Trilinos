// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SCALED_ADJOINT_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_SCALED_ADJOINT_LINEAR_OP_DECL_HPP


#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Concrete decorator <tt>LinearOpBase</tt> subclass that wraps a
<tt>LinearOpBase</tt> object and adds on an extra scaling factor and/or a
transpose enum.

This class represents a scaled, adjointed (transposed) linear operator
<tt>M</tt> of the form:

\verbatim
 
  M = scalar * op(Op)
\endverbatim

where <tt>Op</tt> is another <tt>LinearOpBase</tt> object, <tt>scalar</tt> is
a <tt>Scalar</tt>, and the operation <tt>op(Op)</tt> is specified by a
<tt>EOpTransp</tt> and is given as <tt>op(Op) = Op</tt> (<tt>NOTRANS</tt>), or
<tt>op(Op) = Op^T</tt> (<tt>TRANS</tt>), or <tt>op(Op) = Op^H</tt>
(<tt>CONJTRANS</tt>).  Of course the operator <tt>M</tt> is not constructed
explicitly but instead just applies the decorated operator <tt>Op</tt> by
modifying the <tt>apply()</tt> function that calls <tt>Op.apply()</tt>.

This subclass is designed to allow the efficient handling of multiple implicit
scalings and/or adjoints (transposes) and allow these implicit
transformations to be reversed.  A sequence of scalings/adjoints from some
original <tt>LinearOpBase</tt> object <tt>origOp</tt> is shown as:

\verbatim

 M = scalar_n * op_n( ... scalar_2 * op_2( scalar_1 * op_1( origOp ) ) ... )

   =>

 M  = overallScalar * overall_op(origOp)
\endverbatim

where <tt>overallScalar</tt> and <tt>overall_op(...)</tt> are the cumulative
transformations given as:

\verbatim

 overallScalar = scalar_n * ... * scalar_2 * ... * scalar_1

 overall_op(origOp) = op_n( ... op_2( op_1( ... op_n( origOp ) ... ) ) )
\endverbatim
 
Each individual transformation pair <tt>(scalar_i,op_i(...))</tt> is specified
with arguments <tt>Scalar scalar</tt> and <tt>EOpTransp transp</tt>.  The
overall scaling is returned using <tt>this->overallScalar()</tt>, the overall
adjoint enum is returned using <tt>this->overallTransp()</tt>, and the
original linear operator is returned from <tt>this->getOrigOp()</tt>.

The operator most recently wrapped is returned from <tt>this->getOp()</tt>.
The individual scalings and transformations are not exposed from this
interface but can be viewed by calling <tt>this->description()</tt> with a
verbosity level of ???.  The arguments passed into the constructor
<tt>DefaultScaledAdjointLinearOp()</tt> or <tt>initialize()</tt> can always be
extracted using <tt>uninitialize()</tt>.

This subclass keeps track of all of the individual scalings <tt>scalar_i</tt>
and adjoining operations <tt>op_i(...)</tt> so that the <tt>description()</tt>
function can print this out and also so that the operations can be reversed.

The copy constructor and assignment operators are declared private since some
thought needs to be put into what they should mean.

Note: This class does not maintain any specific cached information about the
original operator so it is safe, in some cases, to manipulate the original
operator and still keep <tt>this</tt> intact and automatically updated for any
changes that are made.

\ingroup Thyra_Op_Vec_ANA_Development_grp

*/
template<class Scalar>
class DefaultScaledAdjointLinearOp
  : virtual public ScaledAdjointLinearOpBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  DefaultScaledAdjointLinearOp();

  /** \brief Calls <tt>initialize()</tt>.
   *
   * Note, instead of calling this constructor directly consider using the
   * non-member functions described belos which create dynamically allocated
   * <tt>RCP</tt>-wrapped objects.
   */
  DefaultScaledAdjointLinearOp(
    const Scalar &scalar,
    const EOpTransp &transp,
    const RCP<LinearOpBase<Scalar> > &Op
    );

  /** \brief Calls <tt>initialize()</tt>.
   *
   * Note, instead of calling this constructor directly consider using the
   * non-member functions described belos which create dynamically allocated
   * <tt>RCP</tt>-wrapped objects.
   */
  DefaultScaledAdjointLinearOp(
    const Scalar &scalar,
    const EOpTransp &transp,
    const RCP<const LinearOpBase<Scalar> > &Op
    );

  /** \brief Initialize with an operator with by defining adjoint (transpose) and
   * scaling arguments.
   *
   * \param scalar [in] Scalar argument defining <tt>scalar_0</tt> (see
   * introduction).
   *
   * \param transp [in] Transpose argument defining <tt>op_0(...)</tt> (see
   * introduction).
   *
   * \param Op [in] Smart pointer to linear operator (persisting
   * relationship).
   *
   *
   * Preconditions:<ul>
   * <li><tt>Op.get() != NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>ToDo: Fill these in!!!!
   * </ul>
   */
  void initialize(
    const Scalar &scalar,
    const EOpTransp &transp,
    const RCP<LinearOpBase<Scalar> > &Op
    );
  
  /** \brief Initialize with an operator with by defining adjoint (transpose) and
   * scaling arguments.
   *
   * \param scalar [in] Scalar argument defining <tt>scalar_0</tt> (see
   * introduction).
   *
   * \param transp [in] Transpose argument defining <tt>op_0(...)</tt> (see
   * introduction).
   *
   * \param Op [in] Smart pointer to linear operator (persisting
   * relationship).
   *
   * Preconditions:<ul>
   * <li><tt>Op.get() != NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>ToDo: Fill these in!!!!
   * </ul>
   */
  void initialize(
 const Scalar &scalar
 ,const EOpTransp &transp
 ,const RCP<const LinearOpBase<Scalar> > &Op
    );

  /** \brief Return the non-const linear operator passed into
   * <tt>initialize()</tt>.
   */
  RCP<LinearOpBase<Scalar> > getNonconstOp();

  /** \brief Return the const linear operator passed into
   * <tt>initialize()</tt>.
   */
  RCP<const LinearOpBase<Scalar> > getOp() const;

  /** \brief Set to uninitialized and (optionally) extract the objects passed into <tt>initialize()</tt>.
   *
   * Postconditions:<ul>
   * <li><tt>this->range().get()==NULL</tt>
   * </ul>
   */
  void uninitialize();

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Outputs
   * <tt>DefaultScaledAdjointLinearOp<Scalar>{this->getOrigOp().description())</tt>
   * along with the dimensions.
   */
  std::string description() const;

  /** \brief Prints out the original operator as well as all of the scalings
   * and transpositions in the order that they occurred.
   *
   * This function outputs different levels of detail based on the value passed in
   * for <tt>verbLevel</tt>:
   *
   * ToDo: Finish documentation!
   */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Return the range space of the logical linear operator.
   *
   * Simply returns: \code

   return ( this->overallTransp()==NOTRANS ? this->getOrigOp()->range() : this->getOrigOp()->domain() );
   \endcode
   */
  RCP<const VectorSpaceBase<Scalar> > range() const;

  /** \brief Return the domain space of the logical linear operator.
   *
   * Simply returns: \code

   return ( this->overallTransp()==NOTRANS ? this->getOrigOp()->domain() : this->getOrigOp()->range() );
   \endcode
   */
  RCP<const VectorSpaceBase<Scalar> > domain() const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** \name Overridden from ScaledAdointLinearOpBase */
  //@{

  /** \brief . */
  Scalar overallScalar() const;
  /** \brief . */
  EOpTransp overallTransp() const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > getNonconstOrigOp();
  /** \brief . */
  RCP<const LinearOpBase<Scalar> > getOrigOp() const;

  //@}

protected:
  
  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Return if the operation is supported on the logical linear
   * operator.
   *
   * Simply returns: \code

   return this->getOrigOp()->opSupported(trans_trans(this->overallTransp(),M_trans));
   \endcode
   */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief Apply the linear operator (or its transpose) to a multi-vector :
   * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
   *
   * Simply calls: \code

   this->getOrigOp()->apply(trans_trans(M_trans,this->overallTransp()),X,Y,(this->overallScalar()*alpha),beta)
   \endcode
   */
 void applyImpl(
   const EOpTransp M_trans,
   const MultiVectorBase<Scalar> &X,
   const Ptr<MultiVectorBase<Scalar> > &Y,
   const Scalar alpha,
   const Scalar beta
   ) const;

  //@}

private:

  // ////////////////////////////////
  // Private types

  template <class Scalar2>
  struct ScalarETransp {
    Scalar2   scalar;
    EOpTransp   transp;
    ScalarETransp()
      {}
    ScalarETransp( const Scalar2 &_scalar, const EOpTransp &_transp )
      : scalar(_scalar), transp(_transp)
      {}
  };

  typedef std::vector<ScalarETransp<Scalar> >  allScalarETransp_t;

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > CNLOC;
  
  // ////////////////////////////////
  // Private data members

  CNLOC origOp_;
  Scalar overallScalar_;
  EOpTransp overallTransp_;
  int my_index_;
  
 RCP<allScalarETransp_t> allScalarETransp_;
  
  // ////////////////////////////////
  // Private member functions

  void initializeImpl(
    const Scalar &scalar,
    const EOpTransp &transp,
    const RCP<const LinearOpBase<Scalar> > &Op,
    const bool isConst
    );
  CNLOC getOpImpl() const;
  void assertInitialized() const;

  // Not defined and not to be called
  DefaultScaledAdjointLinearOp(const DefaultScaledAdjointLinearOp<Scalar>&);
  DefaultScaledAdjointLinearOp<Scalar>& operator=(const DefaultScaledAdjointLinearOp<Scalar>&);

};


/** \brief Build an implicit non-<tt>const</tt> scaled linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstScale(
  const Scalar &scalar,
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit <tt>const</tt> scaled linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
scale(
  const Scalar &scalar,
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit non-<tt>const</tt> adjoined linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstAdjoint(
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit <tt>const</tt> adjoined linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
adjoint(
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit non-<tt>const</tt> transposed linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstTranspose(
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit <tt>const</tt> transposed linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
transpose(
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit non-<tt>const</tt> scaled and/or adjoined
 * (transposed) linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(scale,transp,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstScaleAndAdjoint(
  const Scalar &scalar, const EOpTransp &transp,
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


/** \brief Build an implicit <tt>const</tt> scaled and/or adjoined
 * (transposed) linear operator.
 *
 * Returns <tt>Teuchos::rcp(new
 * DefaultScaledAdjointLinearOp<Scalar>(scale,transp,Op)</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>Op.get()!=NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li><tt>return.get()!=NULL</tt>
 * </ul>
 *
 * \relates DefaultScaledAdjointLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
scaleAndAdjoint(
  const Scalar &scalar, const EOpTransp &transp,
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label = ""
  );


// /////////////////////////////////
// Inline members


template<class Scalar>
inline
DefaultScaledAdjointLinearOp<Scalar>::DefaultScaledAdjointLinearOp()
  :overallScalar_(Teuchos::ScalarTraits<Scalar>::zero())
  ,overallTransp_(NOTRANS)
  ,my_index_(0)
{}


template<class Scalar>
inline
DefaultScaledAdjointLinearOp<Scalar>::DefaultScaledAdjointLinearOp(
  const Scalar &scalar
  ,const EOpTransp &transp
  ,const RCP<LinearOpBase<Scalar> > &Op
  )
  :overallScalar_(Teuchos::ScalarTraits<Scalar>::zero())
  ,overallTransp_(NOTRANS)
{
  this->initialize(scalar,transp,Op);
}


template<class Scalar>
inline
DefaultScaledAdjointLinearOp<Scalar>::DefaultScaledAdjointLinearOp(
  const Scalar &scalar
  ,const EOpTransp &transp
  ,const RCP<const LinearOpBase<Scalar> > &Op
  )
  :overallScalar_(Teuchos::ScalarTraits<Scalar>::zero())
  ,overallTransp_(NOTRANS)
{
  this->initialize(scalar,transp,Op);
}


template<class Scalar>
inline
void DefaultScaledAdjointLinearOp<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( origOp_.getConstObj().get() == NULL );
#endif
}


}	// end namespace Thyra


// /////////////////////////////////
// Inline non-members


template<class Scalar> inline
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstScale(
  const Scalar &scalar,
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        scalar,NOTRANS,Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::scale(
  const Scalar &scalar,
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
    new DefaultScaledAdjointLinearOp<Scalar>(scalar,NOTRANS,Op)
    );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstAdjoint(
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::adjoint(
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        Teuchos::ScalarTraits<Scalar>::one(),CONJTRANS,Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstTranspose(
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::transpose(
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        Teuchos::ScalarTraits<Scalar>::one(),TRANS,Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstScaleAndAdjoint(
  const Scalar &scalar,
  const EOpTransp &transp,
  const RCP<LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(scalar,transp,Op)
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


template<class Scalar> inline
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::scaleAndAdjoint(
  const Scalar &scalar,
  const EOpTransp &transp,
  const RCP<const LinearOpBase<Scalar> > &Op,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> >
    salo = Teuchos::rcp(
      new DefaultScaledAdjointLinearOp<Scalar>(
        scalar, transp, Op
        )
      );
  if (label.length())
    salo->setObjectLabel(label);
  return salo;
}


#endif	// THYRA_DEFAULT_SCALED_ADJOINT_LINEAR_OP_DECL_HPP
