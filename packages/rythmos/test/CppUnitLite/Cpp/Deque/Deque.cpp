#include "Deque.h"

Deque::Deque() : itsSize(0), itsCap(10)
{
	itsData = new int[itsCap];
}

Deque::Deque(const Deque& other) : itsSize(other.size()),
                                   itsCap(other.capacity())
{
	copy_data(other.itsData);
}

int Deque::capacity() const
{
	return itsCap;
}

int Deque::size() const
{
	return itsSize;
}

void Deque::push_back(int i)
{
	if (itsSize >= itsCap)
	{
		grow(0);
	}
	itsData[itsSize++] = i;
}

void Deque::push_front(int i)
{
	itsSize++;
	grow(1);
	itsData[0] = i;
}

int Deque::pop_back()
{
	return itsData[--itsSize];
}

int Deque::pop_front()
{
	int value = itsData[0];
	shrink();
	itsSize--;
	return value;
}

void Deque::grow(int place)
{
	itsCap *= 2;
	array_copy(place);
}

void Deque::array_copy(int place, int start)
{
	int* temp = new int[itsCap];
	array_copy(itsData, temp, start, place);
	delete [] itsData;
	itsData = temp;
}

void Deque::array_copy(const int* source,
					  int* destination,
					  int source_start,
					  int destination_start)
{
	for (int i = 0; i <= itsSize; i++)
		destination[i + destination_start] = source[i + source_start];
}

void Deque::shrink()
{
	array_copy(0, 1);
}

void Deque::copy_data(int* data)
{
	itsData = new int[itsCap];
	array_copy(data, itsData);
}

Deque& Deque::operator=(const Deque& d)
{
	if (this != &d)
	{
		itsCap = d.capacity();
		itsSize = d.size();
		delete [] itsData;
		copy_data(d.itsData);
	}
	return *this;
}