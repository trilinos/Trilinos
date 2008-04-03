#ifndef DEQUE_H
#define DEQUE_H

template <class T>
class Deque
{
	public:
		Deque();
		Deque(const Deque& d);
		int size() const;
		int capacity() const;
		void push_back(const T& i);
		void push_front(const T& i);
		T& pop_back();
		T pop_front();
		Deque& operator=(const Deque& d);

	protected:
		T* itsData;

	private:
        void array_copy(const T* source,
        			   T* destination,
        			   int source_start = 0,
        			   int destination_start = 0);
		void copy_data(int* data);
		void grow(int place);
		void shrink();
		void array_copy(int place, int start = 0);

		int itsSize;
		int itsCap;
};
#include "Deque.cii"
#endif