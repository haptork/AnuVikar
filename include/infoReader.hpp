/*!
 * @file
 * functions to read parcas input file and xyz file and use other
 * functions to get the results and print them.
 * */

#ifndef INFOREADER_ANUVIKAR_HPP
#define INFOREADER_ANUVIKAR_HPP

#include <string>

#include <cluster2features.hpp>
#include <helper.hpp>

namespace avi {

std::pair<std::string, std::string> getInfileFromXyzfile(std::string xyzfile);

std::pair<avi::XyzFileType, bool> getSimulationCode(std::string fname);

std::tuple<avi::InputInfo, avi::ExtraInfo, bool>
extractInfoParcas(std::string fname, std::string ftag);

std::tuple<avi::InputInfo, avi::ExtraInfo, bool>
extractInfoLammps(std::string fname, std::string ftag);

std::tuple<avi::InputInfo, avi::ExtraInfo, bool> infoFromStdIn();

} // namespace avi
#endif