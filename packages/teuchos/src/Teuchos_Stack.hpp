#ifndef TEUCHOS_STACK_H
#define TEUCHOS_STACK_H

#include "Teuchos_Array.hpp"

namespace Teuchos
{
	using std::string;

	/**
	 * \ingroup Containers
	 * Templated LIFO stack.
	 * @author Kevin Long
	 */


	template<class T>
		class Stack
		{
		public:
			/** Construct an empty stack */
			inline Stack();
	
			/** push an element onto the top of the stack */
			inline void push(const T& data);

			/** pop the top element from the top of the stack */
			inline T pop();

			/** peek at the top element */
			inline T peek();

			/** get the number of elements in the stack */
			inline int size() const {return list_.length();}

			/** read elements into an Array */
			inline Array<T> arrayify();

		private:
			Array<T> list_;
		};


	// create an empty list
	template<class T> inline Stack<T>::Stack()
		: list_()
		{;}

	// put a new entry at the beginning of the list

	template<class T> inline void Stack<T>::push(const T& data)
		{
			list_.append(data);
		}

	// return last entry from list, then unhook it from the list.

	template<class T> inline T Stack<T>::pop()
		{
			if (list_.length()==0) Error::raise("Stack<T>::get() on empty Stack");
            T rtn = list_[list_.length()-1];
			list_.remove(list_.length()-1);
			return rtn;
		}

	template<class T> inline T Stack<T>::peek() 
		{
			if (list_.length()==0) Error::raise("Stack<T>::get() on empty Stack");
			return list_[list_.length()-1];
		}


	template<class T> inline Array<T> Stack<T>::arrayify()
		{
          return list_;
		}

	template<class T> inline string toString(const Stack<T>& stack)
		{
			return toString(list_);
		}



}
#endif





