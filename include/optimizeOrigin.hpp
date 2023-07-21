
/*!
 * @file
 * Functions to correct the origin to minimize the overall offset of all the
 * atoms
 * */

#ifndef OPTIMIZEORIGIN_ANUVIKAR_HPP
#define OPTIMIZEORIGIN_ANUVIKAR_HPP

#include <helper.hpp>

namespace avi {

Coords correctOrigin(const std::vector<avi::Coords> &atoms,
                     avi::Coords origin, double latConst);

Coords estimateOrigin(const std::vector<avi::Coords> &atoms,
                      const double &latConst);

} // namespace avi
#endif
