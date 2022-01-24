/*!
 * @file
 * Function to print the results as json
 * */
#ifndef PRINTJSON_ANUVIKAR_HPP
#define PRINTJSON_ANUVIKAR_HPP

#include <array>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <helper.hpp>
#include <results.hpp>

namespace avi {
void printJson(std::ostream &outfile, const InputInfo &i, const ExtraInfo &ei,
               const resultsT &res);

void resToKeyValue(std::ostream &outfile, const resultsT &res);
void infoToKeyValue(std::ostream &outfile, const InputInfo &i,
                    const ExtraInfo &ei);
void configToKeyValue(std::ostream &outfile, const avi::Config &c);
} // namespace avi

#endif