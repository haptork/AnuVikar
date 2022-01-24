/*!
 * @file
 * Finding different features like cluster size distribution, distance and
 * angles of defects from PKA origin.
 * */
#include <cmath>

#include <UnionFind.hpp>
#include <helper.hpp>
#include <results.hpp>

#include <iostream>

// group defects into clusters
avi::DefectVecT avi::groupDefects(const avi::DefectVecT &defects,
                                            const double &latticeConst) {
  // std::cout << "Grouping defects " << defects.size() << '\n' << std::flush;
  using UF = avi::UnionFind<2, avi::DefectT>;
  auto nn = (std::sqrt(3) * latticeConst) / 2 + avi::invars::epsilon;
  auto nn2sqr = (latticeConst * latticeConst) + avi::invars::epsilon; // std::sqrt(3) * info.latticeConst + 0.01;
  auto nn4sqr = (nn * nn) * 4;
  auto pred = [nn2sqr, nn4sqr](const avi::DefectT &a,
                         const avi::DefectT &b) {
    using namespace DefectTWrap;
    if ((isVacancy(a) && isVacancy(b) && (isSurviving(a) && isSurviving(b))) || 
        (isInterstitial(a) && isInterstitial(b) && (isSurviving(a) && isSurviving(b)))) { // both vacancies and surviving
      return calcDistSqr(coords(a), coords(b)) < nn4sqr;
    } 
    return calcDistSqr(coords(a), coords(b)) < nn2sqr;
  };
  UF uf;
  for (const auto &it : defects) {
    uf.uniteIf(it, pred);
  }
  // std::cout << "Finished grouping defects\n" << std::flush;
  return uf.getAll();
  // coords, isInterstitial, ClusterId, isSurviving
}

// cluster id and their sizes
avi::ClusterSizeMapT
avi::clusterSizes(const avi::DefectVecT &defects) {
  avi::ClusterSizeMapT clusterSize;
  using namespace avi::DefectTWrap;
  for (const auto &it : defects) {
    clusterSize[clusterId(it)].all++;
    if (!isSurviving(it)) continue;
    if (isInterstitial(it))
      clusterSize[clusterId(it)].surviving++;
    else
      clusterSize[clusterId(it)].surviving--;
  }
  return clusterSize;
}

avi::ClusterIVMapT
avi::clusterIVType(const avi::ClusterIdMapT &a,
                        avi::ClusterSizeMapT &b) {
  avi::ClusterIVMapT res;
  for (const auto &it : a)
    res[it.first] = b[it.first].surviving;
  return res;
}

// ignore dumbbells or similar defects group from cluster list
void avi::ignoreSmallClusters(avi::DefectVecT &defects,
                                   avi::ClusterSizeMapT &clusterSize) {
  using namespace avi::DefectTWrap;
  using avi::invars::minClusterPoints;
  using avi::invars::minClusterSize;
  for (auto &it : defects) {

    if (abs(clusterSize[clusterId(it)].surviving) < minClusterSize
        && clusterSize[clusterId(it)].all < minClusterPoints) {

      if (clusterSize[clusterId(it)].surviving > 0) {
        clusterId(it, -1);
      } else {
        clusterId(it, 0); // setting clusterId of small ones to zero
      }
    }
  }
}

// cluster ids mapped to defect ids that the clusters have
avi::ClusterIdMapT
avi::clusterMapping(const avi::DefectVecT &defects) {
  using namespace avi::DefectTWrap;
  avi::ClusterIdMapT clusterIds;
  int i = 0;
  for (const auto &it : defects) {
    if (clusterId(it) > 0) {
      clusterIds[clusterId(it)].push_back(
          i); // adding cluster ids and defect index
    }
    ++i;
  }
  return clusterIds;
}

// fraction of defects in cluster
std::tuple<int, double, double>
avi::getNDefectsAndClusterFractions(const avi::DefectVecT &defects) {
  using namespace avi::DefectTWrap;
  auto inClusterI = 0;
  auto inClusterV = 0;
  auto singlesI = 0;
  auto singlesV = 0;
  for (const auto &it : defects) {
    if (!isSurviving(it)) continue;
    if (clusterId(it) <= 0) {
      if (isInterstitial(it))
        singlesI++;
      else
        singlesV++;
    } else {
      if (isInterstitial(it))
        inClusterI++;
      else
        inClusterV++;
    }
  }
  auto nDefects = inClusterI + singlesI;
  // auto nDefects = inClusterV + singlesV;
  double inClusterFractionI =
      (nDefects > 0) ? (double)(inClusterI)*100.0 / nDefects : 0;
  double inClusterFractionV =
      (nDefects > 0) ? (double)(inClusterV)*100.0 / nDefects : 0;
  return std::make_tuple(nDefects, inClusterFractionI, inClusterFractionV);
}

// cluster to cluster features mapping
avi::ClusterFeatMapT
avi::clusterFeatures(const avi::DefectVecT &defects,
                          const avi::ClusterIdMapT &clusters,
                          avi::ClusterSizeMapT &clusterCounts,
                          double latticeConst) {
  using namespace avi::DefectTWrap;
  avi::ClusterFeatMapT clusterFeats;
  using avi::invars::minClusterPoints;
  using avi::invars::minClusterSize;
  for (const auto &it : clusters) {
    if (abs(clusterCounts[it.first].surviving) < minClusterSize && clusterCounts[it.first].all < minClusterPoints) continue;
    std::vector<avi::Coords> clusterCoords;
    std::vector<bool> isI;
    for (const auto &jt : it.second) {
      auto x = coords(defects[jt]);
      clusterCoords.push_back(avi::Coords{{x[0], x[1], x[2]}});
      isI.push_back(isInterstitial(defects[jt]));
    }
    clusterFeats[it.first] =
        avi::pairHists(clusterCoords, isI, latticeConst);
  }
  return clusterFeats;
}

// maximum size of the interstitial and vacancy clusters
std::tuple<int, int>
avi::getMaxClusterSizes(avi::ClusterSizeMapT &clusterCounts,
                             const avi::ClusterIdMapT &clusters) {
  auto maxClusterSizeV = 0;
  auto maxClusterSizeI = 0;
  for (const auto &it : clusters) {
    const int sz = clusterCounts[it.first].surviving;
    if (sz > 0 && sz > maxClusterSizeI)
      maxClusterSizeI = sz;
    else if (sz < 0 && sz < maxClusterSizeV)
      maxClusterSizeV = sz;
  }
  return std::make_tuple(maxClusterSizeI, std::abs(maxClusterSizeV));
}