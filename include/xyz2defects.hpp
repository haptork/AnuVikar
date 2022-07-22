/*!
 * @file
 * implements an O(log(N)) algorithm for detecting defects given a
 * single frame file having xyz coordinates of all the atoms
 * */
#ifndef XYZ2DEFECTS_ANUVIKAR_HPP
#define XYZ2DEFECTS_ANUVIKAR_HPP

#include <string>
#include <tuple>

#include <AddOffset.hpp>
#include <helper.hpp>
#include <xyzReader.hpp>
#include <results.hpp>

namespace avi {


struct dispCoords {
    std::vector<std::pair<Coords, size_t>> sias;
    std::vector<Coords> vacs;
};


// exposed for testing
Coords getInitialMax(const Coords &origin, const Coords &maxes);

Coords getInitialMaxFcc(const Coords &origin, const Coords &maxes);

DefectRes xyz2defectsTime(InputInfo &mainInfo, ExtraInfo &extraInfo,
                      const Config &config, std::istream &infile, frameStatus &fs);

DefectRes displaced2defectsTime(InputInfo &mainInfo, ExtraInfo &extraInfo,
                      const Config &config, std::istream &infile, frameStatus &fs);

DefectRes atoms2defects(std::tuple<xyzFileStatus, std::vector<offsetCoords>, std::vector<std::vector<double>>> atoms,
              InputInfo &info, ExtraInfo &extraInfo, const Config &config, bool isBcc);

DefectRes displacedAtoms2defects(std::tuple<xyzFileStatus, dispCoords,
                                 std::vector<std::vector<double>>> d, double lc);
} // namespace avi

#endif // XYZ2DEFECTS_ANUVIKAR_HPP
