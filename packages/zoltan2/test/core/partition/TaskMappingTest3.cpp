#include "Zoltan2_TaskMapping.hpp"
#include "Zoltan2_TestHelpers.hpp"
#include "Tpetra_MultiVector_decl.hpp"

#include "GeometricGenerator.hpp"
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
        int narg,
        char **argv,
        std::string &procF,
        part_t &nx,
        part_t &ny,
        part_t &nz, bool &divide_prime, int &ranks_per_node, string &taskGraphFile, string &taskCoordFile){

    bool isprocset = false;
    int ispartset = 0;

    for(int i = 0; i < narg; ++i){
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
        else if(identifier == "DP"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, tmp2, '=');
            int val;
            stream >> val;
            if (val) divide_prime = true;
            ispartset++;
        }
        else if(identifier == "RPN"){
          std::stringstream stream(std::stringstream::in | std::stringstream::out);
          stream << tmp;
          getline(stream, tmp2, '=');
          stream >> ranks_per_node;
          ispartset++;
        } else if(identifier == "TG"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, taskGraphFile, '=');

            stream >> taskGraphFile;
        }
        else if(identifier == "TC"){
          std::stringstream stream(std::stringstream::in | std::stringstream::out);
          stream << tmp;
          getline(stream, taskCoordFile, '=');
          stream >> taskCoordFile;
        }

        else {
            throw "Invalid argument at " + tmp;
        }

    }
    if(!(ispartset >= 3&& isprocset)){
        throw "(PROC && PART) are mandatory arguments.";
    }

}
int main(int narg, char *arg[]){

    Tpetra::ScopeGuard tscope(&narg, &arg);

    typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
    typedef inputAdapter_t::part_t part_t;

    //if (narg != 3){
    //    std::cout << "Usage: " << arg[0] << " PART=partGeoParams.txt PROC=procGeoParams.txt" << std::endl;
    //    std::terminate();
    //}
    part_t numParts = 0;
    zscalar_t **partCenters = NULL;
    int coordDim = 0;

    part_t numProcs = 0;
    zscalar_t **procCoordinates = NULL;
    int procDim = 0;


    string taskGraphFile = "";
    string taskCoordFile = "";
    bool divide_prime = false;
    part_t jobX = 1, jobY = 1, jobZ = 1;
    string procfile = "";
    int rank_per_node = 16;

    const RCP<Comm<int> > commN;
    RCP<Comm<int> >comm =  Teuchos::rcp_const_cast<Comm<int> >
            (Teuchos::DefaultComm<int>::getDefaultSerialComm(commN));

    part_t *task_communication_xadj_ = NULL;
    part_t *task_communication_adj_ = NULL;
    zscalar_t *task_communication_adjw_ = NULL;
    try {

        getArgVals<part_t>(
                narg,
                arg,
                procfile ,
                jobX, jobY, jobZ, divide_prime, rank_per_node, taskGraphFile, taskCoordFile);

        coordDim = 3;
        procDim = 3;
        numParts = jobZ*jobY*jobX;

        //numProcs = numParts;
        //std::cout << "part:" << numParts << " proc:" << procfile << std::endl;
        if (taskGraphFile == "" || taskCoordFile == "")
        {

            partCenters = new zscalar_t * [coordDim];
            for(int i = 0; i < coordDim; ++i){
                partCenters[i] = new zscalar_t[numParts];
            }


            task_communication_xadj_ = new part_t [numParts+1];
            task_communication_adj_ = new part_t [numParts * 6];

            int prevNCount = 0;
            task_communication_xadj_[0] = 0;
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
              task_communication_xadj_[i+1] = prevNCount;
            }
        }
        else {
          int ne = 0;
          FILE *f2 = fopen(taskGraphFile.c_str(), "rb");

          fread(&numParts,sizeof(int),1,f2); // write 10 bytes to our buffer
          fread(&ne,sizeof(int),1,f2); // write 10 bytes to our buffer


          std::cout << "numParts:" << numParts << " ne:" << ne << std::endl;

          task_communication_xadj_ = new part_t [numParts+1];
          task_communication_adj_ = new part_t [ne];
          task_communication_adjw_ = new zscalar_t [ne];

          fread((void *)task_communication_xadj_,sizeof(int),numParts + 1,f2); // write 10 bytes to our buffer
          fread((void *)task_communication_adj_,sizeof(int),ne ,f2); // write 10 bytes to our buffer
          fread((void *)task_communication_adjw_,sizeof(double),ne,f2); // write 10 bytes to our buffer
          fclose(f2);

          f2 = fopen(taskCoordFile.c_str(), "rb");
          fread((void *)&numParts,sizeof(int),1,f2); // write 10 bytes to our buffer
          fread((void *)&coordDim,sizeof(int),1,f2); // write 10 bytes to our buffer

          std::cout << "numParts:" << numParts << " coordDim:" << coordDim << std::endl;

          partCenters = new zscalar_t * [coordDim];
          for(int i = 0; i < coordDim; ++i){
              partCenters[i] = new zscalar_t[numParts];
              fread((void *) partCenters[i],sizeof(double),numParts, f2); // write 10 bytes to our buffer
          }
          fclose(f2);
        }




        {
            std::vector < std::vector <zscalar_t> > proc_coords(procDim);
            std::fstream m(procfile.c_str());
            procCoordinates = new zscalar_t * [procDim];

            part_t i = 0;
            zscalar_t a,b, c;
            m >> a >> b >> c;
            while(!m.eof()){
              proc_coords[0].push_back(a);
              proc_coords[1].push_back(b);
              proc_coords[2].push_back(c);
              ++i;
              m >> a >> b >> c;
            }

            m.close();
            numProcs = i;
            for(int ii = 0; ii < procDim; ++ii){
              procCoordinates[ii] = new zscalar_t[numProcs];
              for (int j = 0; j < numProcs; ++j){
                procCoordinates[ii][j] = proc_coords[ii][j];
              }
            }
        }
        std::cout << "numProcs:" << numProcs << std::endl;

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
        Teuchos::RCP<const Teuchos::Comm<int> > tcomm =Tpetra::getDefaultComm();
        part_t *proc_to_task_xadj_ = new part_t[numProcs+1];
        part_t *proc_to_task_adj_ = new part_t[numParts];
/*
        std::cout << "procDim:" << procDim <<
                " numProcs:" << numProcs <<
                " coordDim:" << coordDim <<
                " numParts" << numParts << std::endl;

        for(part_t j = 0; j < numProcs; ++j){
            std::cout << "proc - coord:" << j << " " << procCoordinates[0][j]<< " " << procCoordinates[1][j]<< " " << procCoordinates[2][j] << std::endl;
        }

        for(part_t j = 0; j < numParts; ++j){
            std::cout << "part - coord:" << j << " " << partCenters[0][j]<< " " << partCenters[1][j]<< " " << partCenters[2][j] << std::endl;
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
                task_communication_adjw_,

                proc_to_task_xadj_, /*output*/
                proc_to_task_adj_, /*output*/

                partArraysize,
                Kokkos::View<part_t*, Kokkos::HostSpace,
                  Kokkos::MemoryUnmanaged>(partArray,(partArraysize == -1 ? 0 : partArraysize)),
                machineDimensions, rank_per_node, divide_prime
                );
        
        if (tcomm->getRank() == 0){
            std::cout << "PASS" << std::endl;
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
        delete [] task_communication_adjw_;

        for (int i = 0; i < coordDim; i++) delete [] partCenters[i];
        delete [] partCenters;
        for (int i = 0; i < procDim; i++) delete [] procCoordinates[i];
        delete [] procCoordinates;
    }
    catch(std::string &s){
        std::cerr << s << std::endl;
    }

    catch(char * s){
        std::cerr << s << std::endl;
    }
}

