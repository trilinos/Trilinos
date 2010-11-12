
#ifndef PANZER_BOUNDARYCONDITION_H
#define PANZER_BOUNDARYCONDITION_H

#include <string>
#include <iostream>
#include <functional>

namespace panzer {

  //! Type of boundary condition - determines evaluators to use
  enum BCType {
    BCT_Dirichlet,
    BCT_Neumann
  };

  /*! \brief Stores input information for Dirichlet boundary conditions.

      An object of this type is created for each Dirichlet boundary
      condition parsed in the input file.  A nodal boundary condition
      specifies the solution value at the nodal points of the surface.
      It can be a constant or a complex function of other unknowns at
      the nodes.
  */
  class BoundaryCondition {
    
  public:

    /*! \brief Ctor for a constant value bc.

        Input file example:

	\verbatim
	bc, dirichlet, TEMP, sideset 6, element block 3, -1.0
	\endverbatim	
    */
    BoundaryCondition(int bc_id,
		      BCType bc_type,
		      int sideset_id,
		      int element_block_id,
		      std::string dof_name,
		      double value);

    /*! \brief Ctor that uses a generic strategy.

        Input file example:

	\verbatim
	bc, dirichlet, VELOCITY, sideset 6, element block 3, my strategy
	\endverbatim
	
    */
    BoundaryCondition(int bc_id,
		      BCType bc_type,
		      int sideset_id,
		      int element_block_id,
		      std::string equation_set_name,
		      std::string strategy);

    /*! \brief Constructor that requires a strategy and method function.

        Input file example:

	\verbatim
	bc, nodal, TEMP, sideset 6, element block 3, method function, 1
	bc, dirichlet, TEMP, sideset 6, element block 3, method function, 1
	\endverbatim
    */
    BoundaryCondition(int bc_id,
		      BCType bc_type,
		      int sideset_id,
		      int element_block_id,
		      std::string equation_set_name,
		      std::string strategy,
		      int method_function_id);

    //! Dtor.
    ~BoundaryCondition();

    //! Returns the nodebc id.  This is a unique identifier for this nodebc condition - needed for unique parameter setting in LOCA and for map key comparisons (strict weak ordering).
    int bcID() const;

    //! Returns the boundary condition type (Dirichlet or Neumann).
    BCType bcType() const;

    //! Returns the set id.  Note: this could be a sideset or a nodeset id - call getSetType() to check for the type.
    int sidesetID() const;

    //! Returns the element block id associated with this sideset.
    int elementBlockID() const;

    //! Returns the unknown name/keyword.
    std::string equationSetName() const;

    //! Returns the value if the function is constant.
    double constantValue() const;

    //! Returns the keyword used to construct a variable provider.
    std::string strategy() const;

    //! Returns the method function id used to import extra data.
    int methodFunctionID() const;

    //! Returns true if this is a constant bc.
    bool isConstant() const;

    //! Returns true if this is a strategy bc.
    bool isStrategy() const;

    //! Returns true if this is a method function bc.
    bool isMethodFunction() const;

    //! A unique string identifier for this boundary condition.
    std::string identifier() const;

    //! Print object using an ostream.
    void print(std::ostream& os) const;

  private:

    void setup(int nodebc_id,
	       BCType bc_type,
	       int set_id,
	       int element_block_id,
	       std::string equation_set_name);

  private:

    //! Type of constructor used.
    enum BCFunctionType {
      //! Constant
      BCFT_Constant,
      //! A strategy class
      BCFT_Strategy,
      //! A strategy that requires a method function,
      BCFT_MethodFunction
    };

  private:

    BCFunctionType m_bcf_type;

    int m_bc_id;

    BCType m_bc_type;

    int m_sideset_id;

    int m_element_block_id;

    std::string m_equation_set_name;

    double m_constant_value;

    std::string m_strategy;

    int m_method_function_id;

  };

  std::ostream& 
    operator<<(std::ostream & os, const panzer::BoundaryCondition& bc);

  struct LessBoundaryCondition {
    
    bool operator()(const panzer::BoundaryCondition& left, 
		    const panzer::BoundaryCondition& right) const
    {
      return left.bcID() < right.bcID();
    }
  };

}

#endif
