#include "Zoltan2_TaskMapping.hpp"
#include <Zoltan2_TestHelpers.hpp>
#include "Tpetra_MultiVector_decl.hpp"

#include <GeometricGenerator.hpp>
#include <string>
#include "Teuchos_XMLParameterListHelpers.hpp"
//#include "Teuchos_MPIComm.hpp"
#define partDIM 3
#define procDIM 3
#define nProcs 200;
#define nParts 200;


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

        std::fstream inParam(paramFileName.c_str());
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
    std::istringstream inParam(inp);
    std::string str;
    getline (inParam,str);
    while (!inParam.eof()){
        if(str[0] != param_comment){
            size_t pos = str.find('=');
            if(pos == std::string::npos){
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
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    stream << argumentline;
    getline(stream, argumentid, '=');
    if (stream.eof()){
        return false;
    }
    stream >> argumentValue;
    return true;
}

template <typename part_t>
void getArgVals(
        int argc,
        char **argv,
        std::string &procF,
        part_t &nx,
        part_t &ny,
        part_t &nz){

    bool isprocset = false;
    int ispartset = 0;

    for(int i = 0; i < argc; ++i){
        string tmp = convert_to_string(argv[i]);
        string tmp2 = "";
        string identifier = "";
        double fval = -1;
        if(!getArgumentValue(identifier, fval, tmp)) continue;

        if(identifier == "PROC"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, procF, '=');

            stream >> procF;
            isprocset = true;
        }
        else if(identifier == "NX"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, tmp2, '=');

            stream >> nx;
            ispartset++;
        }
        else if(identifier == "NY"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, tmp2, '=');

            stream >> ny;
            ispartset++;
        }
        else if(identifier == "NZ"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, tmp2, '=');

            stream >> nz;
            ispartset++;
        }

        else {
            throw "Invalid argument at " + tmp;
        }

    }
    if(!(ispartset == 3&& isprocset)){
        throw "(PROC && PART) are mandatory arguments.";
    }

}
int main(int argc, char *argv[]){

    typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
    typedef inputAdapter_t::part_t part_t;

    Teuchos::GlobalMPISession session(&argc, &argv);
    //if (argc != 3){
    //    cout << "Usage: " << argv[0] << " PART=partGeoParams.txt PROC=procGeoParams.txt" << endl;
    //    exit(1);
    //}
    part_t numParts = 0;
    zscalar_t **partCenters = NULL;
    int coordDim = 0;

    part_t numProcs = 0;
    zscalar_t **procCoordinates = NULL;
    int procDim = 0;



    part_t jobX = 1, jobY = 1, jobZ = 1;
    string procfile = "";

    const RCP<Comm<int> > commN;
    RCP<Comm<int> >comm =  Teuchos::rcp_const_cast<Comm<int> >
            (Teuchos::DefaultComm<int>::getDefaultSerialComm(commN));

    part_t *task_communication_xadj_ = NULL;
    part_t *task_communication_adj_ = NULL;
    try {

        getArgVals<part_t>(
                argc,
                argv,
                procfile ,
                jobX, jobY, jobZ);

        coordDim = 3;
        procDim = 3;
        numParts = jobZ*jobY*jobX;
        numProcs = numParts;
        //cout << "part:" << numParts << " proc:" << procfile << endl;
        {
            partCenters = new zscalar_t * [coordDim];
            for(int i = 0; i < coordDim; ++i){
                partCenters[i] = new zscalar_t[numParts];
            }


            task_communication_xadj_ = new part_t [numParts];
            task_communication_adj_ = new part_t [numParts * 6];

            int prevNCount = 0;
            for (part_t i = 0; i < numParts; ++i) {
              int x = i % jobX;
              int y = (i / (jobX)) % jobY;
              int z = (i / (jobX)) / jobY;
              partCenters[0][i] = x;
              partCenters[1][i] = y;
              partCenters[2][i] = z;

              if (x > 0){
              task_communication_adj_[prevNCount++] = i - 1;
              }
              if (x < jobX - 1){
              task_communication_adj_[prevNCount++] = i + 1;
              }
              if (y > 0){
              task_communication_adj_[prevNCount++] = i - jobX;
              }
              if (y < jobY - 1){
              task_communication_adj_[prevNCount++] = i + jobX;
              }
              if (z > 0){
              task_communication_adj_[prevNCount++] = i - jobX * jobY;
              }
              if (z < jobZ - 1){
              task_communication_adj_[prevNCount++] = i + jobX * jobY;
              }
              task_communication_xadj_[i] = prevNCount;
            }
        }



        {
            std::fstream m(procfile.c_str());
            procCoordinates = new zscalar_t * [procDim];
            for(int i = 0; i < procDim; ++i){
                procCoordinates[i] = new zscalar_t[numParts];
            }
            part_t i = 0;
            while(i < numProcs){
                m >> procCoordinates[0][i] >> procCoordinates[1][i] >> procCoordinates[2][i];
                //cout << "i:" <<i << endl;
                ++i;

            }
            m.close();
        }


        /*
        Zoltan2::CoordinateCommunicationModel<zscalar_t,zscalar_t,int> *cm =
                new Zoltan2::CoordinateCommunicationModel<zscalar_t,zscalar_t,int>(
                        procDim, procCoordinates,
                        coordDim, partCenters,
                        numProcs, numParts);

        Zoltan2::Environment *env = new Zoltan2::Environment();
        Zoltan2::CoordinateTaskMapper <inputAdapter_t, int> *ctm=
                new Zoltan2::CoordinateTaskMapper<inputAdapter_t,int>(env, cm);

        */
        RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
        part_t *proc_to_task_xadj_ = new part_t[numProcs];
        part_t *proc_to_task_adj_ = new part_t[numParts];
/*
        cout << "procDim:" << procDim <<
                " numProcs:" << numProcs <<
                " coordDim:" << coordDim <<
                " numParts" << numParts << endl;

        for(part_t j = 0; j < numProcs; ++j){
            cout << "proc - coord:" << j << " " << procCoordinates[0][j]<< " " << procCoordinates[1][j]<< " " << procCoordinates[2][j] << endl;
        }

        for(part_t j = 0; j < numParts; ++j){
            cout << "part - coord:" << j << " " << partCenters[0][j]<< " " << partCenters[1][j]<< " " << partCenters[2][j] << endl;
        }
*/
        /*
        int partArray[3];
        partArray[0] = 8;
        partArray[1] = 4;
        partArray[2] = 16;
        */
        part_t *partArray = NULL;
        int partArraysize = -1;
        //part_t hopper[3];
        //hopper[0] = 17;
        //hopper[1] = 8;
        //hopper[2] = 24;
        part_t *machineDimensions = NULL;
        //machineDimensions = hopper;
        Zoltan2::coordinateTaskMapperInterface<part_t, zscalar_t, zscalar_t>(
                tcomm,
                procDim,
                numProcs,
                procCoordinates,

                coordDim,
                numParts,
                partCenters,

                task_communication_xadj_,
                task_communication_adj_,
		NULL,

                proc_to_task_xadj_, /*output*/
                proc_to_task_adj_, /*output*/
		
                partArraysize,
                partArray,
                machineDimensions
                );

        if (tcomm->getRank() == 0){
            cout << "PASS" << endl;
        }
        /*
        delete ctm;
        delete cm;
        delete env;
        */
        delete [] proc_to_task_xadj_;
        delete [] proc_to_task_adj_;
        delete [] task_communication_xadj_;
        delete [] task_communication_adj_;

        for (int i = 0; i < coordDim; i++) delete [] partCenters[i];
        delete [] partCenters;
        for (int i = 0; i < procDim; i++) delete [] procCoordinates[i];
        delete [] procCoordinates;
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

