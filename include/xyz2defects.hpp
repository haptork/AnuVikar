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

// Gives next expected lattice site given a lattice site
// in this way it enumerates all the lattice sites
// Currently only works for bcc
// TODO: extend for fcc

// exposed for testing
Coords getInitialMax(const Coords &origin, const Coords &maxes);

Coords getInitialMaxFcc(const Coords &origin, const Coords &maxes);

DefectRes xyz2defectsTime(InputInfo &mainInfo, ExtraInfo &extraInfo,
                      const Config &config, std::istream &infile, frameStatus &fs);

DefectRes displaced2defectsTime(InputInfo &mainInfo, ExtraInfo &extraInfo,
                      const Config &config, std::istream &infile, frameStatus &fs);

DefectRes atoms2defects(std::pair<xyzFileStatus, std::vector<offsetCoords>> atoms,
              InputInfo &info, ExtraInfo &extraInfo, const Config &config);

DefectRes atoms2defectsFcc(std::pair<xyzFileStatus, std::vector<offsetCoords>> atoms,
              InputInfo &info, ExtraInfo &extraInfo, const Config &config);

DefectRes displacedAtoms2defects(std::pair<xyzFileStatus, std::array<std::vector<Coords>, 2>> d,
                           double lc);
} // namespace avi

#endif // XYZ2DEFECTS_ANUVIKAR_HPP
