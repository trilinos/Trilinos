#ifndef OPTIONSFORTESTING_HPP_
#define OPTIONSFORTESTING_HPP_

#include <string>

#include <stk_unit_test_utils/getOption.h>

struct Options
{
    std::string mMeshFileName;
    bool mOverRideTest;
    bool mDebugZoltan;
    std::string mZoltanDebugLevel;

    int mNumSubdomains;
    std::string mParmetisMethod;
    std::string mParmetisOutputLevel;
    std::string mParmetisIter;
    std::string mParmetisCheckGraph;
    std::string mLargestDegreeFirst;
    std::string mOutputFilename;

    bool mDeleteFiles;

    double mToleranceForFaceSearch;
    double mToleranceForParticleSearch;

    int mNumTargetProcs;

    std::string getMeshFileName() const
    {
        return mMeshFileName;
    }

    bool overRideTest() const
    {
        return mOverRideTest;
    }

    bool debugZoltan() const
    {
        return mDebugZoltan;
    }

    std::string getZoltanDebugLevel() const
    {
        return mZoltanDebugLevel;
    }

    int numSubdomains() const
    {
        return mNumSubdomains;
    }

    std::string getStringNumSubdomains() const
    {
        std::ostringstream os;
        os << mNumSubdomains;
        return os.str();
    }

    std::string getPartmetisMethod() const
    {
        return mParmetisMethod;
    }

    std::string getPartmetisOutputLevel() const
    {
        return mParmetisOutputLevel;
    }

    std::string getParmetisIter() const
    {
        return mParmetisIter;
    }

    std::string getParmetisCheckGraph() const
    {
        return mParmetisCheckGraph;
    }

    std::string getLargestDegreeFirst() const
    {
        return mLargestDegreeFirst;
    }

    bool deleteFiles() const
    {
        return mDeleteFiles;
    }

    double getToleranceForFaceSearch() const
    {
        return mToleranceForFaceSearch;
    }

    double getToleranceForParticleSearch() const
    {
        return mToleranceForParticleSearch;
    }

    int getNumTargetProcs() const
    {
        return mNumTargetProcs;
    }

    std::string getOutputFilename() const
    {
        return mOutputFilename;
    }

    //////////////////////////////////////////////

    void setNumTargetProcs(int num_target_procs)
    {
        mNumTargetProcs = num_target_procs;
    }

    void setMeshFileName(const std::string& filename)
    {
        mMeshFileName = filename;
    }

    void setOverRideTest(bool override)
    {
        mOverRideTest = override;
    }

    void setDebugZoltan(bool debug)
    {
        mDebugZoltan = debug;
        if(mDebugZoltan == true)
        {
            this->setZoltanDebugLevel("1");
        }
    }

    void setZoltanDebugLevel(const std::string& debuglevel)
    {
        mZoltanDebugLevel = debuglevel;
    }

    void setNumSubdomains(int numsub)
    {
        mNumSubdomains = numsub;
    }

    void setPartmetisMethod(const std::string& method)
    {
        mParmetisMethod = method;
    }

    void setPartmetisOutputLevel(const std::string& outputlevel)
    {
        mParmetisOutputLevel = outputlevel;
    }

    void setParmetisIter(const std::string& iterlevel)
    {
        mParmetisIter = iterlevel;
    }

    void setParmetisCheckGraph(const std::string& checkgraph)
    {
        mParmetisCheckGraph = checkgraph;
    }

    void setLargestDegreeFirst(const std::string& degreefirst)
    {
        mLargestDegreeFirst = degreefirst;
    }

    void setDeleteFiles(bool deleteFiles)
    {
        mDeleteFiles = deleteFiles;
    }

    void setToleranceForFaceSearch(double tol)
    {
        mToleranceForFaceSearch = tol;
    }

    void setToleranceForParticleSearch(double tol)
    {
        mToleranceForParticleSearch = tol;
    }

    void setOutputFilename(const std::string& filename)
    {
        mOutputFilename=filename;
    }
    //////////////////////////////////

    Options() :
            mMeshFileName(), mOverRideTest(false), mDebugZoltan(false),
                    mZoltanDebugLevel("0"),
                    mNumSubdomains(1), mParmetisMethod("PartKway"),
                    mParmetisOutputLevel("0"), mParmetisIter("100"),
                    mLargestDegreeFirst("L"), mOutputFilename("subdomain.exo"), mDeleteFiles(true),
                    mToleranceForFaceSearch(0.1), mToleranceForParticleSearch(1.0),
                    mNumTargetProcs(0)
    {
    }
};

inline Options getOptionsForTest(const std::string &defaultMesh)
{
    Options local;
    const std::string generatedMeshSpec = unitTestUtils::getOption("-i", defaultMesh);
    local.setMeshFileName(generatedMeshSpec);

    std::string manualRun = unitTestUtils::getOption("-manual", "no");
    if(manualRun != "no")
    {
        local.setOverRideTest(true);
    }

    std::string debugZoltanLevel = unitTestUtils::getOption("-zdl", "0");
    local.setZoltanDebugLevel(debugZoltanLevel);

    std::string debugZoltan = unitTestUtils::getOption("-z", "no");
    if(debugZoltan != "no")
    {
        local.setDebugZoltan(true);
    }

    {
        std::string targetProcs = unitTestUtils::getOption("-t", "0");
        int numTargetProcs = 0;
        std::istringstream is(targetProcs);
        is >> numTargetProcs;
        local.setNumTargetProcs(numTargetProcs);
    }

    const std::string outputFilename = unitTestUtils::getOption("-o", "subdomain.exo");
    local.setOutputFilename(outputFilename);

    {
        const std::string numDecompProcs = unitTestUtils::getOption("-n", "3");
        int numsub = 1;
        std::istringstream is(numDecompProcs);
        is >> numsub;
        local.setNumSubdomains(numsub);
    }

    std::string parmetisMethod = unitTestUtils::getOption("-m", "PartKway");
    local.setPartmetisMethod(parmetisMethod);

    std::string parmetisOutputLevel = unitTestUtils::getOption("-o", "0");
    local.setPartmetisOutputLevel(parmetisOutputLevel);

    std::string parmetisIter = unitTestUtils::getOption("-iter", "100");
    local.setParmetisIter(parmetisIter);

    std::string parmetisCheckGraph = unitTestUtils::getOption("-check", "1");
    local.setParmetisCheckGraph(parmetisCheckGraph);

    std::string largestDegreeFirst = unitTestUtils::getOption("-v", "L");
    local.setLargestDegreeFirst(largestDegreeFirst);

    std::string deleteFile = unitTestUtils::getOption("-d", "yes");
    if(deleteFile != "yes")
    {
        local.setDeleteFiles(false);
    }

    {
        std::string toleranceForFaceSearch = unitTestUtils::getOption("-tolFace", "0.1");
        std::istringstream stream(toleranceForFaceSearch);
        double faceSearchTolerance = 0;
        stream >> faceSearchTolerance;
        local.setToleranceForFaceSearch(faceSearchTolerance);
    }

    {
        std::string toleranceForParticleSearch = unitTestUtils::getOption("-tolPart", "1.0");
        std::istringstream stream(toleranceForParticleSearch);
        double particleSearchTolerance = 0;
        stream >> particleSearchTolerance;
        local.setToleranceForParticleSearch(particleSearchTolerance);
    }

    return local;
}



#endif /* OPTIONSFORTESTING_HPP_ */
