  // Scenario: We are solving a system of equations Ax = b.
  //           We know A and b, and have an approximation to x
  //           These two versions of code attempt compute the 
  //           2-norm of the residual r where r = b - Ax.



  // Code fragment 1:  THIS CODE WILL NOT WORK
  // Fails because the statement A.Multiply() typicall involves 
  // interprocessor communication but processor 0 is the only one executing
  // the Multiply() method.  The user program will probably stall.
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      A.Multiply(false, q, r); // Compute Ax, store in r
      r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 2:  THIS CODE WILL NOT WORK
  // The Multiply() method will work.  The Update() method will complete, but results 
  // will be incorrect.  Finally this code will stall on the call to Norm2() 
  // because it involved a collective operation where all processor must 
  // participate.
  A.Multiply(false, q, r); // Compute Ax, store in r
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 3:  THIS CODE WILL NOT WORK
  // The Multiply() and Update() methods will work.  
  // As with the previous segment this code will stall on the call to Norm2() 
  // because it involved a collective operation where all processor must 
  // participate.
  A.Multiply(false, q, r); // Compute Ax, store in r
  r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 4:  THIS CODE WILL WORK
  // All methods will work.  Only the output statement itself is restricted to
  // processor 0.
  A.Multiply(false, q, r); // Compute Ax, store in r
  r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
  if (Comm().MyPID()==0)  // Only processor 0 will execute this code
	 cout << "2-norm of b - Ax = " << r << endl;
