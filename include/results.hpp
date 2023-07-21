/*!
 * @file
 * Finding different results like cluster size distribution, distance and
 * angles of defects from PKA origin.
 * */
#ifndef RESULTS_ANUVIKAR_HPP
#define RESULTS_ANUVIKAR_HPP

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <cluster2features.hpp>
#include <helper.hpp>

namespace avi {

// coords, isInterstitial, ClusterId, isSurviving
using DefectT = std::tuple<avi::Coords, bool, int, bool>;

using DefectVecT =
    std::vector<std::tuple<avi::Coords, bool, int, bool>>;

using DefectRes = std::tuple<xyzFileStatus, ErrorStatus, DefectVecT, std::vector<int>, std::vector<std::vector<double>>>;

struct ClusterSizeT {
  int surviving;
  int all;
};

using ClusterSizeMapT = std::unordered_map<int, ClusterSizeT>;

using ClusterIdMapT = std::unordered_map<int, std::vector<int>>;

using ClusterIVMapT = std::unordered_map<int, int>;

namespace DefectTWrap {
static inline auto coords(const DefectT &_d) { return std::get<0>(_d); }
static inline auto isInterstitial(const DefectT &_d) { return std::get<1>(_d); }
static inline auto clusterId(const DefectT &_d) { return std::get<2>(_d); }
static inline auto isSurviving(const DefectT &_d) { return std::get<3>(_d); }
static inline void coords(DefectT &_d, avi::Coords c) {
  std::get<0>(_d) = c;
}
static inline void isInterstitial(DefectT &_d, bool a) { std::get<1>(_d) = a; }
static inline void clusterId(DefectT &_d, int id) { std::get<2>(_d) = id; }
static inline void isSurviving(DefectT &_d, bool s) { std::get<3>(_d) = s; }
static inline auto isVacancy(const DefectT &_d) { return !std::get<1>(_d); }
} // namespace DefectTWrap

struct resultsT {
  ErrorStatus err;
  int nDefects;
  int nSia;
  int nVac;
  int nClusters;
  int maxClusterSizeI;
  int maxClusterSizeV;
  double inClusterFractionI;
  double inClusterFractionV;
  DefectVecT defects;
  std::vector<int> coDefects;
  std::vector<std::vector<double>> extraCols;
  std::array<std::vector<double>, 2> dists;
  std::array<std::vector<double>, 2> angles;
  std::unordered_map<int, std::vector<int>> clusters;
  std::unordered_map<int, int> clustersIV;
  std::unordered_map<int, featT> feats;
};

// group defects into clusters
DefectVecT groupDefects(const DefectVecT &defects, const double &latticeConst, const avi::Coords &box);

// surviving and all defects count for groups of all sizes, including single
// defects that are currently given a clusterId of their own
ClusterSizeMapT clusterSizes(const DefectVecT &defects);

// change clusterId to zero for defects that belong to clusters having lesser
// defects than the input given criterion
void ignoreSmallClusters(DefectVecT &defects, ClusterSizeMapT &);

// for each cluster a list of coord indices of defects that belong to it
ClusterIdMapT clusterMapping(const DefectVecT &defects);

// surviving defects count for clusters ignoring single defects labelled as 0
// clusterId
ClusterIVMapT clusterIVType(const ClusterIdMapT &, ClusterSizeMapT &);

// number of defects and fraction of interstitials and vacancies in cluster
std::tuple<int, int, double, double>
getNDefectsAndClusterFractions(const avi::DefectVecT &defects);

using ClusterFeatMapT = std::unordered_map<int, featT>;

// cluster features based on pair distances, triad angles and neighbourhood
// count
ClusterFeatMapT clusterFeatures(const avi::DefectVecT &defects,
                                const avi::ClusterIdMapT &clusters,
                                avi::ClusterSizeMapT &clusterCounts,
                                double latticeConst, avi::Coords box);

// maximum of surviving defects among all of clusters
// (interstitial and vacancy cluster separately)
std::tuple<int, int>
getMaxClusterSizes(avi::ClusterSizeMapT &clusterCounts,
                   const avi::ClusterIdMapT &clusters);

} // namespace avi
#endif // ANUVIKAR_RESULTS_HPP
