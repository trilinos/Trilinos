#include "IOContFileUtils.H"

#include "iostream"
#include "sstream"

using namespace std;

bool WriteHeaderToContFile( const string & fileName,
    const Teuchos::ParameterList & fileParams )
{
  // The file to open
  ofstream oFile(fileName.c_str());

  // Writing the header
  oFile << "#" << setw(6) << "ID";

  // Looping on the parameters
  Teuchos::map<string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i) 
    oFile << setw(15) << fileParams.name(i);
  oFile << std::endl;

  // Closing
  oFile.close();

  return true;
}

bool UpdateContFile( const string & fileName,
    const int & idStep,
    const Teuchos::ParameterList & fileParams )
{

  // The file to open
  ofstream oFile(fileName.c_str(), ios_base::app);

  // Writing the id
  oFile << scientific << setw(7)  << idStep;

  // Looping on the parameters
  Teuchos::map<string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i) 
    oFile << scientific << setw(15) << fileParams.entry(i);
  oFile << std::endl;

  // Closing
  oFile.close();

  return true;
}

bool RestartContFile( const string & fileName, const int & idStep )
{
  // Copying the continuation in a backup file
  string fileNameBak = fileName + ".bak";
  string command = "cp " + fileName + " " + fileNameBak;
  system (command.c_str());

  // String of the line to cut from
  ostringstream ostream;
  ostream << idStep + 1; 
  string lineNumber = ostream.str();

  // Cutting the file
  command =  "sed '" + lineNumber + ",$ d' " + fileNameBak + " > " + fileName;
  system (command.c_str());

  return true;
}

bool TouchContFileParameters( Teuchos::ParameterList & fileParams )
{
  // Either int or double type
  int dummyInt;
  double dummyDouble;

  // Looping the list
  Teuchos::map<string, Teuchos::ParameterEntry>::const_iterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i) {

    if (fileParams.isType<int>(fileParams.name(i)))
      dummyInt = fileParams.get<int>(fileParams.name(i));

    if (fileParams.isType<double>(fileParams.name(i)))
      dummyDouble = fileParams.get<double>(fileParams.name(i));
  }

  return true;
}
