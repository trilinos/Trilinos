#include"Isorropia_EpetraMatcher.hpp"
using namespace std;

int main(int argc, char** argv) {

	vector<int> mt;
	
	if(argc>1)
	{	
		matcher pm(argv[1]);
		mt=pm.get_matching();
		/*for(int i=0;i<32;i++)
			cout<<mt[i]<<" ";
		cout<<endl;*/
	}
	else
		cout<<"Specify input file.."<<endl;
	
	return 0;
}
