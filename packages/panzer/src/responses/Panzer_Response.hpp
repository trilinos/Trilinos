#ifndef __Panzer_Response_hpp__
#define __Panzer_Response_hpp__

#include "Panzer_LinearObjFactory.hpp"

namespace panzer {

/** A simple structure describing the type 
  * value of the response.
  */
struct ResponseId {
   std::string name;
   std::string type;

   ResponseId(const std::string & n,
              const std::string & t)
      : name(n), type(t) {}

   std::string getString() const 
   { return name + "("+type+")"; }

   //! (*this) == id
   bool operator==(const ResponseId & id) const
   { return id.name==name && id.type==type; }

   //! (*this) < id
   bool operator<(const ResponseId & id) const
   { 
      if(name==id.name)
         return type<id.type;
      return name<id.name;
   }
};

ResponseId buildResponse(const std::string & n,
                         const std::string & t)
{ return ResponseId(n,t); }

/** Provide access to responses through a user specified
  */
template <typename TraitsT>
class Response {
public:
   Response(const ResponseId & rid) 
      : rid_(rid), value_(0.0), derivative_(Teuchos::null)
      , hasValue_(false), hasDerivative_(false)
   {}

   Response(const Response & r) 
      : rid_(r.rid_), value_(r.value_), derivative_(r.derivative_)
      , hasValue_(r.hasValue), hasDerivative_(r.hasDerviative_)
   {}
   
   //! What are the details of this response.
   const ResponseId & getResponse() const
   { return rid_; }

   //! Get the value for this response
   typename TraitsT::RealType getValue() const
   { return value_; }
 
   //! set value for 
   void setValue(typename TraitsT::RealType v)
   { value_ = v; hasValue_ = true; }

   /** The derivative is assumed to be stored in the residual vector.
     *
     * \returns A pointer to the linear object container containing
     *          the desired residual object. If no derivative is available
     *          <code>Teuchos::null</code> is returned.
     * \note One could also (and likely should) use Thyra here
     */
   Teuchos::RCP<LinearObjContainer> getDerivative() const
   { return derivative_; }

   //! Set the derivative container for this response
   void setDerivative(const Teuchos::RCP<LinearObjContainer> & loc)
   { derivative_ = loc; hasDerivative_ = true;}

   bool hasValue() const { return hasValue_; }
   bool hasDerivative() const { return hasDerivative_; }

private:
   Response(); // hide me

   ResponseId rid_;
   typename TraitsT::RealType value_;
   Teuchos::RCP<LinearObjContainer> derivative_;
   bool hasValue_, hasDerivative_;
};

}

#endif
