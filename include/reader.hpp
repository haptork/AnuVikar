/*!
 * @file
 * functions to read parcas input file and xyz file and use other
 * functions to get the results and print them.
 * */

#ifndef PARCASREADER_ANUVIKAR_HPP
#define PARCASREADER_ANUVIKAR_HPP

#include <string>

#include <xyzReader.hpp>
#include <cluster2features.hpp>
#include <helper.hpp>
#include <results.hpp>
#include <infoReader.hpp>

namespace avi {


std::pair<avi::xyzFileStatus, avi::ErrorStatus> processTimeFile(avi::InputInfo &info,
                                     avi::ExtraInfo &extraInfo,
                                     const avi::Config &config, std::istream &infile, avi::frameStatus &fs, std::ostream &outfile, bool isFirst);

std::pair<avi::ErrorStatus,int> processFileTimeCmd(std::string xyzfileName,
                                            std::ostream &outfile,
                                            const Config &config, int id, avi::InputInfo, avi::ExtraInfo, avi::infoFrom);

}
#endif