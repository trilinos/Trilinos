#ifndef UNITTESTUTILS_HPP
#define UNITTESTUTILS_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <Geom_AxisAlignedBB.h>
#include <Geom_Search.h>
#include <search/ContactRangeSearch.h>
#include <search/ContactCommunication.h>

#include "gtkTraitsForSearch.h"
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/OctTreeOps.hpp>

typedef stk::search::IdentProc<int,int> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> StkBox;
typedef std::pair<StkBox,Ident> StkBoxWithId;
typedef std::vector< StkBoxWithId > StkBoxVector;

typedef geometry::AxisAlignedBB GtkBox;
typedef std::vector<std::pair<Ident,Ident> > SearchResults;
typedef std::pair<GtkBox,Ident> BoxWithId;
typedef std::vector< BoxWithId > GtkBoxVector;

template<class VolumeType>
VolumeType generateBoundingVolume(double x, double y, double z, double radius);

template<>
inline
Point generateBoundingVolume<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
inline
Sphere generateBoundingVolume<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x,y,z),radius);
}

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
inline
StkBox generateBoundingVolume< StkBox >(double x, double y, double z, double radius)
{
  Point min_corner(x-radius,y-radius,z-radius);
  Point max_corner(x+radius,y+radius,z+radius);
  return StkBox(min_corner,max_corner);
}

template <typename VolumeType>
std::pair<VolumeType, Ident> generateBoundingVolume(double x, double y, double z, double radius, int id, int proc)
{
  return std::make_pair(generateBoundingVolume<VolumeType>(x,y,z,radius), Ident(id,proc));
}

inline void gtk_search(GtkBoxVector& local_domain, GtkBoxVector& local_range, MPI_Comm comm, SearchResults& searchResults)
{
    int num_procs = -1;
    int proc_id   = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);

    std::vector<geometry::AxisAlignedBB> rangeBoxes(local_range.size());
    std::vector<geometry::AxisAlignedBB> domainBoxes(local_domain.size());

    for (size_t i=0;i<local_domain.size();i++)
    {
        domainBoxes[i] = local_domain[i].first;
    }

    for (size_t i=0;i<local_range.size();i++)
    {
        rangeBoxes[i] = local_range[i].first;
    }

    std::vector<int> ghost_indices;
    std::vector<int> ghost_procs;
    ACME::BoxA_BoxB_Ghost(domainBoxes, rangeBoxes, comm, ghost_indices, ghost_procs);

    std::vector< std::vector<geometry::AxisAlignedBB> > send_list(num_procs);
    std::vector< std::vector<geometry::AxisAlignedBB> > recv_list(num_procs);

    // i am sending proc 'ghost proc[i]' my range box 'ghost_indices[i]'
    // ghost_indices.size() is total number of communications that need to occur with all procs

    std::vector< std::vector<int> > send_indices(num_procs);
    std::vector< std::vector<int> > recv_indices(num_procs);

    for (size_t i=0;i<ghost_indices.size();i++)
    {
        send_list[ghost_procs[i]].push_back(rangeBoxes[ghost_indices[i]]);
        int id = local_range[ghost_indices[i]].second.id();
        send_indices[ghost_procs[i]].push_back(id);
    }

    ACME::Parallel_Data_Exchange(send_indices, recv_indices, comm );
    ACME::Parallel_Data_Exchange(send_list, recv_list, comm);

    for (size_t i=0;i<recv_list.size();i++)
    {
        for (size_t j=0;j<recv_list[i].size();j++)
        {
            rangeBoxes.push_back(recv_list[i][j]);
            local_range.push_back(std::make_pair(recv_list[i][j], Ident(recv_indices[i][j], i)));
        }
    }

    std::vector<int> interaction_list;
    std::vector<int> first_interaction;
    std::vector<int> last_interaction;

    geometry::BoxA_BoxB_Search(domainBoxes, rangeBoxes, interaction_list, first_interaction, last_interaction);

    typedef std::set <std::pair<Ident,Ident> > localJunk;
    localJunk localResults;

    // Ident box1, Ident box2
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Ident box1_ident = local_domain[i].second;
        for (int j=first_interaction[i];j<last_interaction[i];j++)
        {
            //            Ident box2_ident = local_range[j].second;
            Ident box2_ident = local_range[interaction_list[j]].second;
            localResults.insert(std::make_pair(box1_ident, box2_ident));
        }
    }

    localJunk tmp;
    stk::search::communicate< std::pair<Ident,Ident>, std::pair<Ident,Ident> >(comm, localResults, tmp);
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(searchResults));
//    std::sort(searchResults.begin(), searchResults.end());
}

extern int gl_argc;
extern char** gl_argv;

inline std::string getOption(const std::string& option, const std::string defaultString="no")
{
    std::string returnValue = defaultString;
    if ( gl_argv != 0 )
    {
        for (int i=0;i<gl_argc;i++)
        {
            std::string input_argv(gl_argv[i]);
            if ( option == input_argv )
            {
                if ( (i+1) < gl_argc )
                {
                    returnValue = std::string(gl_argv[i+1]);
                }
                break;
            }
        }
    }
    return returnValue;
}

#include <exodusMeshInterface.h>

inline void createBoundingBoxForElement(const sierra::Mesh::LocalNodeId *connectivity, const int numNodesPerElement,
        const std::vector<double> &coordinates, std::vector<double>& boxCoordinates)
{
    int spatialDim = 3;
    double *minCoordinates = &boxCoordinates[0];
    double *maxCoordinates = &boxCoordinates[spatialDim];

    int firstNode=0;
    for (int j=0;j<spatialDim;j++)
    {
        minCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
        maxCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
    }

    for (int i=1;i<numNodesPerElement;i++)
    {
        sierra::Mesh::LocalNodeId nodeId = connectivity[i];
        for (int j=0;j<spatialDim;j++)
        {
           minCoordinates[j] = std::min(minCoordinates[j], coordinates[spatialDim*nodeId+j]);
           maxCoordinates[j] = std::max(maxCoordinates[j], coordinates[spatialDim*nodeId+j]);
        }
    }
    bool inflateBox = true;
    double percentInflation = 10;
    if ( inflateBox )
    {
        for (int i=0;i<spatialDim;i++)
        {
            double dist = maxCoordinates[i]-minCoordinates[i];
            if ( dist <= 1e-8 ) dist = 0.001;
            double inflation = dist*(0.5*percentInflation)/100.0;
            minCoordinates[i] -= inflation;
            maxCoordinates[i] += inflation;
        }
    }
}
#endif
