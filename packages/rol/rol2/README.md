ROL 2.0 Design Changes
----------------------

  - New components in namespace ROL2

  - Parallel directory tree rol/rol2/src

  - Functions that returned `std::vector<std::string>` outputs are now type void and
    write to streams

  - Separation of classes in declarations and definitions 
    ROL2_ClassName_Decl.hpp and ROL2_ClassName_Def.hpp except for abstract base 
    classes or functions that do not rely on complete types.

  - Filenames for class in namespaces include the namespace in their name, e.g.
    ROL2_TypeU_Algorithm_Decl.hpp

  - Algorithms and Steps are no longer separate entities. Former ROL::Step code 
    refactored into Algorithm subtypes

  - Four general classes of Algorithm, which inherit from ROL2::Algorithm
    and live in namespaces TypeU (Unconstrained), TypeB (Bound), TypeE (Equality), 
    and TypeG (General). 

  - Use of member types for data containers associated with an object
    - AlgorithmState is now Algorithm::State
    - SecantState is now Secant::State
  
  - The above data containers are private members of the class that 
    primarily uses them. Access is granted following this standard design:

    class Algorithm {
    public:
      struct State { \* ... *  };
      const State& getState() const { return *state_; }
      void setState( const Ptr<State>& state ) { state_ = state; }
    protected:
      State& getState() { return *state_; }
    private:
      Ptr<State> state_;
    };

    The getState methods allow us to use covariant return types with read access
    granted to other objects/functions and write access to derived types without unnecessary
    changes to the reference count. The setState allows the state_ member to be 
    set to a type derived from Algorithm::State if needed. 
    
  - Use of scoped enums as member types
    - TRFlag is replaced by TrustRegion::Flag
      SUCCESS is replaced by TrustRegion::Flag::Success

  - string-enum conversion functions become static member functions,
    for example TrustRegion::flagToString

  - utility functions that are specific to a class type become static member functions
    of that class

  - Use of C++14 variable templates for type-dependent constants instead of functions.
    For example ROL_EPSILON<Real>() is replaced by ROL_EPSILON<Real>;

  Minor things:
  - Deletion of void as a function argument
  - Deletion of const modifier on objects that are passed by value
  - Consistent code formatting
    - Default argument values set in function declaration
    - Long function signatures vertically align arguments  
  - use of auto when type is clear from LHS of assignment
  - class ordering: public first, then protected, then private, except as needed
    for type dependence
  - member variable ordering: Decreasing value of sizeof when possible unless it
    compromises clarity
  - consistent use of override 4
  - Always pass Ptr by reference for multi-threading efficiency
  
    
