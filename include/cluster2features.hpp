/*!
 * @file
 * Function to find the histogram of angles, distances and adjacency as the
 * characteristic feature of a cluster shape
 * */
#ifndef CLUSTER2FEATURES_ANUVIKAR_HPP
#define CLUSTER2FEATURES_ANUVIKAR_HPP

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <helper.hpp>

namespace avi {
// using defectsT = std::vector<std::tuple<Coords, bool, int, bool>>;
constexpr size_t distBins = 40;
constexpr size_t angleBins = 36;
constexpr double angleBinSize = 5.0;
constexpr double maxAngle = 180.0;
constexpr size_t adjBins = (size_t)(maxAngle / angleBinSize);

using distFeatT = std::array<double, distBins>;
using angleFeatT = std::array<double, angleBins>;
using adjNn2FeatT = std::array<double, adjBins>;
using featT = std::tuple<distFeatT, angleFeatT, adjNn2FeatT>;
featT pairHists(const std::vector<avi::Coords> &v,
                const std::vector<bool> &v2, double latConst, avi::Coords box);
} // namespace avi
#endif // CLUSTER2FEATURES_ANUVIKAR_HPP
