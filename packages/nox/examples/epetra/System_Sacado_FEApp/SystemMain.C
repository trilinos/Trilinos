#include <iostream>

#include "AppModelEval.H"

int main(int argc, char *argv[]) {

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    /* Construct first model  g1(p1) */
    AppModelEval App1(30);
    Teuchos::RCP<Epetra_Vector> p1 = Teuchos::rcp((new Epetra_Vector(*(App1.get_p_init(0)))));
    Teuchos::RCP<Epetra_Vector> g1 = Teuchos::rcp((new Epetra_Vector(*(App1.get_g_map(0)))));
    EpetraExt::ModelEvaluator::InArgs params_in1 = App1.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs responses_out1 = App1.createOutArgs();


    /* Construct second model  g2(p2) */
    AppModelEval App2(20);
    Teuchos::RCP<Epetra_Vector> p2 = Teuchos::rcp((new Epetra_Vector(*(App2.get_p_init(0)))));
    Teuchos::RCP<Epetra_Vector> g2 = Teuchos::rcp((new Epetra_Vector(*(App2.get_g_map(0)))));
    EpetraExt::ModelEvaluator::InArgs params_in2 = App2.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs responses_out2 = App2.createOutArgs();


    /* Evaluate first model */
    params_in1.set_p(0,p1);
    responses_out1.set_g(0,g1);
    App1.evalModel(params_in1, responses_out1);

    cout << "Finished eval of first model: Params, Responses " << endl;
    p1->Print(cout); g1->Print(cout);


    /* Couple response of first model to input of second model */
    (*p2)[1] = (*g1)[0];

    /* Evaluate second model */
    params_in2.set_p(0,p2);
    responses_out2.set_g(0,g2);
    App2.evalModel(params_in2, responses_out2);

    cout << "Finished eval of second model: Params, Responses " << endl;
    p2->Print(cout); g2->Print(cout);

#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

  return 0;
}
