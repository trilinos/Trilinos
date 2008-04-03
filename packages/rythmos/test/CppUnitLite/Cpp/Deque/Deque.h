#ifndef DEQUE_H
#define DEQUE_H

class Deque
{
	public:
		Deque();
		Deque(const Deque& d);
		int size() const;
		int capacity() const;
		void push_back(int i);
		void push_front(int i);
		int pop_back();
		int pop_front();
		Deque& operator=(const Deque& d);

	protected:
		int* itsData;

	private:
        void array_copy(const int* source,
        			   int* destination,
        			   int source_start = 0,
        			   int destination_start = 0);
		void copy_data(int* data);
		void grow(int place);
		void shrink();
		void array_copy(int place, int start = 0);

		int itsSize;
		int itsCap;
};
#endif