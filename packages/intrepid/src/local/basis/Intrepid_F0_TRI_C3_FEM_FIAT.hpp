#ifndef INTREPID_F0_TRI_C3_FEM_FIAT_HPP
#define INTREPID_F0_TRI_C3_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F0_TRI_C3_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 3 for 0-forms on TRI cells. 
  Reconstruction space type is COMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     0      |     0      |     0      |     1      |               L_0(u)=u(0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     0      |     1      |     0      |     1      |               L_1(u)=u(1.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     0      |     2      |     0      |     1      |               L_2(u)=u(0.0,1.0)  |
  |========|============|============|============|============|==================================|
  |    3   |     1      |     0      |     0      |     2      |    L_3(u)=u(0.333333333333,0.0)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     0      |     1      |     2      |    L_4(u)=u(0.666666666667,0.0)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     1      |     0      |     2      |L_5(u)=u(0.666666666667,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    6   |     1      |     1      |     1      |     2      |L_6(u)=u(0.333333333333,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    7   |     1      |     2      |     0      |     2      |    L_7(u)=u(0.0,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    8   |     1      |     2      |     1      |     2      |    L_8(u)=u(0.0,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    9   |     2      |     0      |     0      |     1      |L_9(u)=u(0.333333333333,0.333333333333)  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_0                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TRI                   |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_COMPLETE    |
  |-----------------------|-----------------------------------|
  |  degree               |  3                          |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F0_TRI_C3_FEM_FIAT: public Basis<Scalar> {
  private:
  
  /** \brief Dimension of the space spanned by the basis = number of degrees of freedom.
  */
  int numDof_;
  
  /**\brief Lookup table for the DoF's local enumeration (DoF Id) by its local DoF tag
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > tagToEnum_;
  
  /**\brief Lookup table for the DoF's local DoF tag by its local enumeration (DoF Id)
  */
  Teuchos::Array<LocalDofTag> enumToTag_;
  
  /**\brief "true" if both lookup arrays have been set by initialize()
  */
  bool isSet_;

  /**\brief coefficients of nodal basis in terms of orthogonal polynomials */
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats1_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm_data() {
    static Scalar vdm_data[] = { 3.3333333333333312e-02 , 3.3333333333333354e-02 , 3.3333333333333347e-02 , 7.4999999999999956e-02 , 7.4999999999999956e-02 , 7.5000000000000011e-02 , 7.5000000000000039e-02 , 7.5000000000000011e-02 , 7.4999999999999900e-02 , 4.4999999999999996e-01 , -4.9999999999999933e-02 , 5.0000000000000100e-02 , 1.4274296030894868e-17 , -4.5000000000000007e-01 , 4.4999999999999996e-01 , 4.4999999999999996e-01 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , -4.4999999999999990e-01 , 3.5685740077237171e-17 , -1.6666666666666680e-02 , -1.6666666666666663e-02 , 3.3333333333333326e-02 , -1.4999999999999994e-01 , -1.5000000000000011e-01 , -1.5000000000000011e-01 , 3.0000000000000004e-01 , 2.9999999999999999e-01 , -1.4999999999999999e-01 , 0.0000000000000000e+00 , 2.1428571428571430e-01 , 2.1428571428571425e-01 , -8.4898136324008704e-34 , -2.1428571428571416e-01 , -2.1428571428571438e-01 , 3.2142857142857140e-01 , -7.6469443022651034e-18 , 7.6469443022651065e-18 , 3.2142857142857140e-01 , -6.4285714285714279e-01 , 1.2857142857142856e-01 , -1.2857142857142856e-01 , 4.9960036108132039e-17 , 1.9285714285714295e-01 , -1.9285714285714287e-01 , 1.2857142857142853e-01 , 3.2142857142857140e-01 , -3.2142857142857151e-01 , -1.2857142857142850e-01 , 0.0000000000000000e+00 , 4.2857142857142858e-02 , 4.2857142857142844e-02 , 1.2857142857142856e-01 , 1.1785714285714283e-01 , 1.1785714285714288e-01 , -9.6428571428571475e-02 , 3.2142857142857133e-02 , 3.2142857142857126e-02 , -9.6428571428571433e-02 , -3.2142857142857134e-01 , -2.2500000000000009e-01 , 2.2499999999999992e-01 , 2.6742912942062758e-33 , 6.7500000000000004e-01 , -6.7499999999999993e-01 , -2.4087874552135090e-17 , 2.4087874552135090e-17 , -2.4087874552135100e-17 , 2.4087874552135090e-17 , -2.6742912942062758e-33 , -1.6071428571428570e-01 , -1.6071428571428575e-01 , 6.3673602243006575e-34 , 1.6071428571428575e-01 , 1.6071428571428570e-01 , 3.2142857142857140e-01 , 5.7352082266988318e-18 , -5.7352082266988341e-18 , 3.2142857142857140e-01 , -6.4285714285714279e-01 , -9.6428571428571405e-02 , 9.6428571428571433e-02 , 3.5685740077237171e-17 , -3.2142857142857077e-02 , 3.2142857142857105e-02 , -3.2142857142857140e-01 , 3.2142857142857140e-01 , -3.2142857142857151e-01 , 3.2142857142857140e-01 , -3.5685740077237171e-17 , -3.2142857142857133e-02 , -3.2142857142857140e-02 , 1.2857142857142859e-01 , -3.2142857142857140e-02 , -3.2142857142857119e-02 , 1.2857142857142864e-01 , -1.9285714285714289e-01 , -1.9285714285714284e-01 , 1.2857142857142856e-01 , 1.2857142857142853e-01 };
    return vdm_data;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 2.0000000000000027e+00 , 2.6010939434075110e-15 , 1.3322676295501878e-15 , -3.3578582091724130e-16 , -8.8817841970012523e-16 , 4.4408920985006262e-16 , 3.3321826744581155e-31 , 9.0449190087829040e-16 , 3.8064789415719794e-16 , 4.4408920985006262e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 2.8033131371785187e-16 , 6.0000000000000018e+00 , -2.2870594307278226e-15 , 2.8081278798873539e-15 , -1.1102230246251565e-16 , 8.6914602499226518e-16 , 9.8278528172711195e-16 , -6.6470072076881765e-15 , 1.3322676295501886e-15 , -1.7953892341081102e-15 , 1.3333333333333348e+00 , -4.4408920985005907e-16 , 3.3333333333333339e+00 , 1.1600697726695504e-16 , -8.8817841970012523e-16 , -8.8817841970012523e-16 , -1.1102230246251484e-16 , -6.9604186360173046e-16 , -2.6645352591003718e-15 , -2.6645352591003757e-15 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 9.9999999999999978e-01 , 4.9880734463515956e-15 , -7.9999999999999893e-01 , 9.9999999999999982e+00 , -3.4972025275692431e-15 , 2.0000000000000023e-01 , 3.3544595672602937e-15 , 2.4397717406456506e-15 , 7.3750529492956830e-15 , 8.3266726846886741e-16 , -1.2064423534260041e-15 , 2.4000000000000017e+00 , -2.2056430755886452e-15 , 2.2463606938302111e-15 , 8.3999999999999986e+00 , 8.1839297243797199e-16 , 4.7854577443575090e-16 , -8.1798230320692499e-15 , 2.7891867720314334e-30 , -2.5122761014374968e-15 , 1.0000000000000002e+00 , -2.7755575615628874e-15 , 3.1999999999999984e+00 , 1.8851133805880218e-16 , 1.1102230246251565e-16 , 4.1999999999999993e+00 , -4.4408920985006193e-16 , -1.7972018431279050e-15 , -4.2188474935755901e-15 , -2.6645352591003757e-15 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 1.0000000000000013e+00 , 1.3005469717037555e-15 , 6.6613381477509392e-16 , -1.6789291045862065e-16 , -4.4408920985006262e-16 , 2.2204460492503131e-16 , 1.6660913372290578e-31 , 4.5224595043914520e-16 , 1.9032394707859897e-16 , 2.2204460492503131e-16 , 3.0000000000000044e+00 , 3.9016409151112663e-15 , 1.7763568394002505e-15 , -6.7021218506963541e-16 , -1.3322676295501878e-15 , 8.8817841970012523e-16 , -1.7763568394002503e-15 , 2.4470221767247951e-17 , -7.6129578831439114e-16 , -8.8817841970012523e-16 , 6.6666666666666785e-01 , 3.0000000000000018e+00 , -3.3333333333333365e-01 , 7.2524186203000977e-16 , 1.1102230246251565e-16 , 5.5511151231257827e-16 , -3.9678577883656926e-16 , -2.8032848151625862e-15 , 9.5161973539299188e-16 , -5.5511151231257827e-17 , 6.6666666666666741e-01 , 5.0000000000000027e+00 , 1.6666666666666643e+00 , 2.5091323576687884e-15 , -6.6613381477509392e-16 , 0.0000000000000000e+00 , -6.9190684927531499e-17 , -6.0352233414910328e-15 , -8.8817841970012208e-16 , -2.6645352591003757e-15 , -1.3333333333333359e+00 , -6.5344555163652028e-15 , 6.6666666666666634e+00 , 1.9027863185310339e-15 , 1.3322676295501878e-15 , -2.6645352591003757e-15 , 5.5989000268007571e-31 , -2.3128891092599169e-15 , -4.3140094671148889e-15 , -4.4408920985006262e-15 , 4.9999999999999989e-01 , 2.4000000000000030e+00 , -4.0000000000000024e-01 , 5.0000000000000009e+00 , -6.0000000000000175e-01 , 1.0000000000000037e-01 , 3.1744646593929406e-15 , 7.7937939548845317e-16 , 4.4646825918854511e-15 , -5.5511151231257827e-17 , 4.9999999999999922e-01 , 1.2000000000000046e+00 , 7.9999999999999960e-01 , 6.9999999999999982e+00 , 4.1999999999999957e+00 , -6.9999999999999984e-01 , 4.2305246607453118e-15 , 2.8436720098849136e-16 , 5.2656292025078871e-15 , -1.5543122344752192e-15 , 4.9999999999999867e-01 , -3.6000000000000019e+00 , 1.5999999999999996e+00 , 1.1139048866966120e-16 , 8.3999999999999986e+00 , 2.1000000000000001e+00 , -5.0423950729136106e-16 , -2.2410927988664288e-15 , -1.7763568394002461e-15 , -8.8817841970012523e-16 , 2.5000000000000013e+00 , 9.5161973539301259e-17 , -2.0000000000000000e+00 , -1.3809814969572098e-15 , 4.4408920985006262e-16 , 1.0500000000000004e+01 , 2.6645352591003761e-15 , -5.9589521525798981e-16 , -3.6795963101862295e-15 , 8.8817841970012523e-16 };
    return dmats1_data;
  }
#endif


  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C3_FEM_FIAT() : numDof_(10), isSet_(false) , vdm_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_vdm_data(),10,10,10) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_dmats0_data(),10,10,10) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_dmats1_data(),10,10,10) ) ) { }
  
  /** \brief Initializes arrays needed for the lookup of the local enumeration (DoF Id) of a 
    degree-of-freedom by its local DoF tag and the reverse lookup of the DoF tag by DoF Id.
   */
  void initialize();
  
  
  /** \brief Returns FieldContainer with multi-indexed quantity with the values of an operator 
  applied to a set of FEM basis functions, evaluated at a set of
  points in a <strong>reference</strong> cell. The rank of the return
  FieldContainer argument depends on the number of points, number of
  basis, functions, their rank, and the rank of the operator applied to them; see 
  FieldContainer::resize(int,int,EField,EOperator,int) for summary of
  the admissible index combinations.
  In particular, the admissible range of the <var>operatorType</var> argument depends on the rank of 
  the basis functions and the space dimension.Because FEM
  reconstruction relies on COMPLETE or INCOMPLETE polynomials, i.e.,
  smooth local spaces, the admissible range of <var>operatorType</var>
  always includes VALUE, D0, D1,...,D10. If derivative order exceeds
  polynomial degree, output container is filled with 0. 
 
      \param outputValues   [out]         - FieldContainer with the computed values (see
                                            implementation for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param operatorType    [in]         - the operator being applied to the basis function
  */    

  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const EOperator                        operatorType) const;
  
  
  /** \brief This method is intended for FVD reconstructions and should not be used here. Its 
    invocation will throw an exception. 
  */  
  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const Cell<Scalar>&                    cell) const;

  
  int getNumLocalDof() const;
  
  
  int getLocalDofEnumeration(const LocalDofTag dofTag);

  
  LocalDofTag getLocalDofTag(int id);

  
  const Teuchos::Array<LocalDofTag> & getAllLocalDofTags();

  
  ECell getCellType() const;

  
  EBasis getBasisType() const;

  
  ECoordinates getCoordinateSystem() const;

  
  int getDegree() const;

};

} // namespace Intrepid

#include "Intrepid_F0_TRI_C3_FEM_FIATDef.hpp"

#endif

