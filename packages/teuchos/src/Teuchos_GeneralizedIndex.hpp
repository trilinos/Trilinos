#ifndef TEUCHOS_GENERALIZEDINDEX_H
#define TEUCHOS_GENERALIZEDINDEX_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Error.hpp"
#include "Teuchos_Stack.hpp"

namespace Teuchos
{
  using std::exception;
  using std::string;

  /**
   *
   */
  class GeneralizedIndex
    {
    public:
      /** */
      GeneralizedIndex(const GeneralizedIndex& index, int i);

      /** */
      GeneralizedIndex(int i);

      /** */
      GeneralizedIndex();

      /** */
      int index() const {return indices_.peek();}

      /** */
      GeneralizedIndex remainder() const ;

    private:
      /** */
      GeneralizedIndex(const Stack<int>& stack);


      /** */
      mutable Stack<int> indices_;



    };
}
#endif
