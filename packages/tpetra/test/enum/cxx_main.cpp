
#include <iostream>

namespace Tpetra {
	enum DataAccess {copy, view};
	
	
	void getArg(DataAccess cv) {
		cout << "Argument = " << cv << endl;
	}
	
	int main(int argc, char* argv[]) {
		cout << "Starting..." << endl;
		int copy = 4;
		getArg(copy);
		return(0);
	}

}
