#include "Zoltan2_TaskMapping.hpp"
#include <Zoltan2_TestHelpers.hpp>
#include "Tpetra_MultiVector_decl.hpp"
#include "Zoltan2_MachineRepresentation.hpp"
#include <GeometricGenerator.hpp>
#include <string>
#include "Teuchos_XMLParameterListHelpers.hpp"
#define partDIM 3
#define procDIM 3
#define nProcs 200;
#define nParts 200;

template <typename scalar_t, typename procId_t>
void randomFill(scalar_t **&partCenters, procId_t &numCoords, int &coorDim){
    partCenters = Zoltan2::allocMemory<scalar_t *>(coorDim);
    for(int i = 0; i < coorDim; ++i){
        partCenters[i] = Zoltan2::allocMemory<scalar_t>(numCoords);
        for(procId_t j = 0; j < numCoords; ++j){
            partCenters[i][j] = i;
        }
    }
}

template <typename scalar_t, typename procId_t>
void getPartCenters(scalar_t **&partCenters, procId_t &numCoords, int &coorDim){
    numCoords = nParts;
    coorDim = partDIM;
    randomFill(partCenters, numCoords, coorDim);
}

template <typename scalar_t, typename procId_t>
void getProcCenters(scalar_t **&procCenters, procId_t &numCoords, int &coorDim){
    numCoords = nProcs;
    coorDim = procDIM;
    randomFill(procCenters, numCoords, coorDim);
}



string trim_right_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

string trim_left_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( s.find_first_not_of( delimiters ) );
}

string trim_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

const char param_comment = '#';
void readGeoGenParams(string paramFileName, Teuchos::ParameterList &geoparams, const RCP<const Teuchos::Comm<int> > & comm){
    std::string input = "";
    char inp[25000];
    for(int i = 0; i < 25000; ++i){
        inp[i] = 0;
    }

    bool fail = false;
    if(comm->getRank() == 0){

        fstream inParam(paramFileName.c_str());
        if (inParam.fail())
        {
            fail = true;
        }
        if(!fail)
        {
            std::string tmp = "";
            getline (inParam,tmp);
            while (!inParam.eof()){
                if(tmp != ""){
                    tmp = trim_copy(tmp);
                    if(tmp != ""){
                        input += tmp + "\n";
                    }
                }
                getline (inParam,tmp);
            }
            inParam.close();
            for (size_t i = 0; i < input.size(); ++i){
                inp[i] = input[i];
            }
        }
    }



    int size = input.size();
    if(fail){
        size = -1;
    }
    comm->broadcast(0, sizeof(int), (char*) &size);
    if(size == -1){
        throw "File " + paramFileName + " cannot be opened.";
    }
    comm->broadcast(0, size, inp);
    istringstream inParam(inp);
    string str;
    getline (inParam,str);
    while (!inParam.eof()){
        if(str[0] != param_comment){
            size_t pos = str.find('=');
            if(pos == string::npos){
                throw  "Invalid Line:" + str  + " in parameter file";
            }
            string paramname = trim_copy(str.substr(0,pos));
            string paramvalue = trim_copy(str.substr(pos + 1));
            geoparams.set(paramname, paramvalue);
        }
        getline (inParam,str);
    }
}
string convert_to_string(char *args){
    string tmp = "";
    for(int i = 0; args[i] != 0; i++)
        tmp += args[i];
    return tmp;
}
bool getArgumentValue(string &argumentid, double &argumentValue, string argumentline){
    stringstream stream(stringstream::in | stringstream::out);
    stream << argumentline;
    getline(stream, argumentid, '=');
    if (stream.eof()){
        return false;
    }
    stream >> argumentValue;
    return true;
}

void getArgVals(
        int argc,
        char **argv,
        std::string &procF,
        std::string &partF){


    bool ispartset = false;

    for(int i = 0; i < argc; ++i){
        string tmp = convert_to_string(argv[i]);
        string identifier = "";
        long long int value = -1; double fval = -1;
        if(!getArgumentValue(identifier, fval, tmp)) continue;
        value = (long long int) (fval);

        if(identifier == "PART"){
            stringstream stream(stringstream::in | stringstream::out);
            stream << tmp;
            getline(stream, partF, '=');

            stream >> partF;
            ispartset = true;
        }

        else {
            throw "Invalid argument at " + tmp;
        }

    }
    if(!(ispartset)){
        throw "(PART) is mandatory argument.";
    }

}
int main(int argc, char *argv[]){
    Teuchos::GlobalMPISession session(&argc, &argv);
    if (argc != 2){
        cout << "Usage: " << argv[0] << " PART=partGeoParams.txt " << endl;
        exit(1);
    }
    zoltan2_partId_t numParts = 0;
    scalar_t **partCenters = NULL;
    int coordDim = 0;

    //zoltan2_partId_t numProcs = 0;
    //scalar_t **procCoordinates = NULL;
    //int procDim = 0;



    string partfile = "";
    string procfile = "";

    const RCP<Comm<int> > commN;
    RCP<Comm<int> >serialcomm =  Teuchos::rcp_const_cast<Comm<int> >
            (Teuchos::DefaultComm<int>::getDefaultSerialComm(commN));
    RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();


    try {

        getArgVals(
                argc,
                argv,
                procfile ,
                partfile);
        //cout << "part:" << partfile << " proc:" << procfile << endl;
        {
            Teuchos::ParameterList geoparams;
            //getPartCenters(partCenters, numParts, coordDim);
            readGeoGenParams(partfile, geoparams, serialcomm);
            GeometricGen::GeometricGenerator<scalar_t, lno_t, gno_t, node_t> *gg = new GeometricGen::GeometricGenerator<scalar_t, lno_t, gno_t, node_t>(geoparams,serialcomm);
            coordDim = gg->getCoordinateDimension();
            numParts = gg->getNumLocalCoords();
            partCenters = new scalar_t * [coordDim];
            for(int i = 0; i < coordDim; ++i){
                partCenters[i] = new scalar_t[numParts];
            }
            gg->getLocalCoordinatesCopy(partCenters);
            /*
        for(int i = 0; i < coordDim; ++i){
            for(int j = 0; j < numParts; ++j){
                cout << partCenters[i][j] << " ";
            }
            cout << endl;
        }
             */
            delete gg;
        }

        const Zoltan2::MachineRepresentation<scalar_t> machine(tcomm);

        typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
        typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;


        const Zoltan2::Environment *env = new Zoltan2::Environment();
        Zoltan2::CoordinateTaskMapper <inputAdapter_t, zoltan2_partId_t> *ctm=
                new Zoltan2::CoordinateTaskMapper<inputAdapter_t,zoltan2_partId_t>(
                        env, tcomm.getRawPtr(), &machine,
                        coordDim, numParts, partCenters
                        );

        ctm->writeMapping2(tcomm->getRank());
        if (tcomm->getRank() == 0){
            cout << "PASS" << endl;
        }
        delete ctm;
        delete env;
    }
    catch(std::string &s){
        cerr << s << endl;
    }

    catch(char * s){
        cerr << s << endl;
    }
    catch(char const * s){
        cerr << s << endl;
    }
}

