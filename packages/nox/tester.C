#include <iostream>
#include "NLS_ParameterList.H"

int main() {
  NLS_ParameterList b;
  NLS_ParameterList& a = b.sublist("Options List");

  b.setParameter("Search Method", "LineSearch");
  a.setParameter("Maximum Inner Iterations", 10);
  a.setParameter("Convergence Tolerance", static_cast<double>(0));

  //  cout << "Search Method A:" << a.getParameter("Search Method", "TrustRegion") << endl;
  //  cout << "Search Method A:" << a.getParameter("Maximum Inner Iterations", 5) << endl;
  //  cout << "Search Method B:" << b.getParameter("Search Method", "TrustRegion") << endl;

  cout << "\nList A:" << endl;
  a.print(cout);
  cout << "\nList B:" << endl;
  b.print(cout);
  //cout << "\nList C:" << endl;
  //c.print(cout);

}

