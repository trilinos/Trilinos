#include <iostream>

void f( int& i )
{
	i = 5;
	std::cout << "\ni = "<<i<<std::endl;
}

void g( int* i )
{
	f(*i);
}

int main()
{
	int a  = 2, *b = NULL;
	g(&a);
	g(b);
	return 0;
}
