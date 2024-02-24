from PyROL.PyROL import ROL

paramlist = ROL.ParameterList
paramlist.sublist("Status Test").set("Gradient Tolerance",1e-8)
paramlist.sublist("Status Test").set("Constraint Tolerance",1e-8)
paramlist.sublist("Status Test").set("Step Tolerance",1e-12)
paramlist.sublist("Status Test").set("Iteration Limit", 250)
paramlist.sublist("Step").set("Type","Line Search")
paramlist.sublist("General").set("Output Level",iprint)