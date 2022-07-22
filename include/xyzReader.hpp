/*!
 * @file
 * functions to read parcas input file and xyz file and use other
 * functions to get the results and print them.
 * */

#ifndef XYZREADER_ANUVIKAR_HPP
#define XYZREADER_ANUVIKAR_HPP

#include <map>
#include <string>

#include <helper.hpp>

namespace avi {

constexpr auto maxColumnsTry = 50;

enum class frameStatus : bool { prelude, inFrame };

enum class lineStatus : int { garbage, coords, inFrameCoords, frameBorder };

std::tuple<avi::lineStatus, avi::Coords, std::vector<double>>
getCoord(const std::string &line, const avi::frameStatus &fs,
         const avi::InputInfo &info, const avi::ExtraInfo &ei);

std::tuple<avi::lineStatus, avi::Coords, std::vector<double>>
getCoordLammps(const std::string &line, const avi::frameStatus &fs, int colStart, int ecStart = -1, int ecEnd = -1);

std::tuple<avi::lineStatus, avi::Coords, std::vector<double>>
getCoordParcas(const std::string &line, const avi::frameStatus &fs, int colStart, int ecStart = -1, int ecEnd = -1);

std::tuple<avi::lineStatus, std::array<avi::Coords, 2>, std::vector<double>>
getCoordDisplaced(const std::string &line, int ecStart = -1, int ecEnd = -1);

std::tuple<avi::lineStatus, avi::Coords, std::vector<double>>
getCoordGeneric(const std::string &line, const avi::frameStatus &fs, int columnStart, int ecStart, int ecEnd);

std::tuple<avi::lineStatus, avi::Coords, std::vector<double>>
getCoordCdb(const std::string &line, const avi::frameStatus &fs,
             const std::string &substrate, int columnStart, int ecStart, int ecEnd);

} // namespace avi
#endif