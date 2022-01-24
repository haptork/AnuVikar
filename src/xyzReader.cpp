/*!
 * @file
 * functions to read parcas input file and xyz file and use other
 * functions to get the results and print them.
 * */

#include <string>

#include <cluster2features.hpp>
#include <helper.hpp>
#include <xyzReader.hpp>

std::pair<avi::lineStatus, avi::Coords>
getCoordGeneric(const std::string &line, const avi::frameStatus &fs, int columnStart) {
  avi::Coords c;
  auto first = std::begin(line);
  auto second = std::begin(line);
  for (auto i = 1; i < columnStart; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) {
      return std::make_pair(avi::lineStatus::frameBorder, c);
    }
  }
  const auto maxTry = columnStart > 0 ? 3 : avi::maxColumnsTry;  // if column start is given then three index after that are coordinates else try max
  for (auto i = 0; i < maxTry; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) {
      if (i > 2) return std::make_pair(avi::lineStatus::coords, c);
      return std::make_pair(avi::lineStatus::frameBorder, c);
    }
    try {
      auto curVal = std::stod(std::string{first, second});
      if (i > 2) {
        c[0] = c[1]; c[1] = c[2]; c[2] = curVal; 
      } else {
        c[i] = curVal;
      }
    } catch (const std::invalid_argument &) {
      return std::make_pair(avi::lineStatus::frameBorder, c);
    } catch (const std::out_of_range &) {
      return std::make_pair(avi::lineStatus::frameBorder, c);
    }
  }
  //std::cout << "from c: ";
  //std::cout << c[0] << ", " << c[1] << ", " << c[2] << '\n';
  return std::make_pair(avi::lineStatus::coords, c);
}

std::pair<avi::lineStatus, avi::Coords>
getCoordCdb(const std::string &line, const avi::frameStatus &fs,
             const std::string &substrate, int columnStart) {
  using avi::lineStatus;
  avi::Coords c;
  auto first = std::find_if(begin(line), end(line),
                            [](int ch) { return !std::isspace(ch); });
  if (first == std::end(line))
    return std::make_pair(lineStatus::garbage, c); // possibly blank line
  auto second =
      std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  std::string word{first, second};
  if (word.size() != substrate.size()) {
    // we can return garbage here and continue with adding in the same frame
    // however frameBorder is assumed for anything else to handle the files with
    // multiple frames.
    return std::make_pair(lineStatus::frameBorder, c);
  }
  for (auto i = 0; i < word.size(); i++) {
    if (std::tolower(word[i]) != std::tolower(substrate[i])) {
      // we can return garbage here ...
      return std::make_pair(lineStatus::frameBorder, c);
    }
  }
  if (columnStart <= 1) columnStart = 2;
  for (int i = 1; i < columnStart; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) return std::make_pair(lineStatus::garbage, c);
  }
  for (int i = 0; i < 3 && first < second; ++i) {
    try {
      c[i] = std::stod(std::string{first, second});
    } catch (const std::invalid_argument &) {
      return std::make_pair(lineStatus::garbage, c);
    } catch (const std::out_of_range &) {
      return std::make_pair(lineStatus::garbage, c);
    }
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second && i < 2) return std::make_pair(lineStatus::garbage, c);

  }
  return std::make_pair(lineStatus::inFrameCoords, c);
}

std::pair<avi::lineStatus, avi::Coords>
avi::getCoord(const std::string &line, const avi::frameStatus &fs,
                   const avi::InputInfo &info,
                   const avi::ExtraInfo &extraInfo) {
  return (info.xyzFileType == avi::XyzFileType::cascadesDbLikeCols)
             ? getCoordCdb(line, fs, extraInfo.substrate, info.xyzColumnStart)
             : (info.xyzFileType == avi::XyzFileType::lammpsWithStdHeader)
                   ? getCoordLammps(line, fs, info.xyzColumnStart)
                   : (info.xyzFileType == avi::XyzFileType::parcasWithStdHeader)
                      ? getCoordParcas(line, fs, info.xyzColumnStart)
                      : getCoordGeneric(line, fs, info.xyzColumnStart);
}

std::pair<avi::lineStatus, avi::Coords>
avi::getCoordLammps(const std::string &line,
                         const avi::frameStatus &fs, int columnStart) {
  avi::Coords c;
  auto first = std::find_if(begin(line), end(line),
                            [](int ch) { return !std::isspace(ch); });
  if (first == std::end(line))
    return std::make_pair(avi::lineStatus::garbage,
                          c); // possibly blank line
  auto second =
      std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  std::string word{first, second};
  if (word == "ITEM:" || word == "TIMESTEP:") {
    return std::make_pair(avi::lineStatus::frameBorder, c);
  }
  if (fs != avi::frameStatus::inFrame) { // Not checking always in favor of
                                              // efficiency but if the xyz file
                                              // has more than a single frame it
                                              // will be a mess
    auto first_temp = first;
    auto second_temp = second;
    while (first_temp != std::end(line) && first_temp < second_temp) {
      std::string word{first_temp, second_temp};
      if (word == "ITEM:" || word == "TIMESTEP:")
        return std::make_pair(avi::lineStatus::frameBorder, c);
      first_temp = std::find_if(second_temp, end(line),
                                [](int ch) { return !std::isspace(ch); });
      second_temp = std::find_if(first_temp, end(line),
                                 [](int ch) { return std::isspace(ch); });
    }
  }
  for (auto i = 1; i < columnStart; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) {
      return std::make_pair(avi::lineStatus::garbage, c);
    }
  }
  auto maxTry = columnStart > 0 ? 3 : avi::maxColumnsTry;  // if column start is given then three index after that are coordinates else try max
  auto curColumn = 0;
  for (auto i = 0; i < maxTry && first < second; ++i) {
    try {
      auto curVal = std::stod(std::string{first, second});
      if (curColumn > 2) {
        c[0] = c[1]; c[1] = c[2]; c[2] = curVal; 
      } else {
        c[curColumn] = curVal;
      }
      curColumn++;
    } catch (const std::invalid_argument &) {
      curColumn = 0;
    } catch (const std::out_of_range &) {
      curColumn = 0;
    }
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  }
  if (curColumn > 2) return std::make_pair(avi::lineStatus::coords, c);
  return std::make_pair(avi::lineStatus::garbage, c);
}

std::pair<avi::lineStatus, avi::Coords>
avi::getCoordParcas(const std::string &line,
                         const avi::frameStatus &fs, int columnStart) {
  avi::Coords c;
  auto first = std::find_if(begin(line), end(line),
                            [](int ch) { return !std::isspace(ch); });
  if (first == std::end(line))
    return std::make_pair(avi::lineStatus::garbage,
                          c); // possibly blank line
  if (std::isdigit(*first) && columnStart != 1)
    return std::make_pair(avi::lineStatus::garbage,
                          c); // possibly first line with frame number
  auto second =
      std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  std::string word{first, second};
  if (word == "Frame" || word == "boxsize") {
    return std::make_pair(avi::lineStatus::frameBorder, c);
  }
  if (columnStart == 0) columnStart = 2;
  for (auto i = 2; i < columnStart; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) {
      return std::make_pair(avi::lineStatus::garbage, c);
    }
  }
  int start = 0;
  if (columnStart == 1) {
    try {
      c[0] = std::stod(word);
    } catch (const std::invalid_argument &) {
      return std::make_pair(avi::lineStatus::garbage, c);
    } catch (const std::out_of_range &) {
      return std::make_pair(avi::lineStatus::garbage, c);
    }
    start = 1;
  }
  auto maxTry = 3;  // three index after column start
  for (int i = start; i < maxTry; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) {
      if (i > 2) return std::make_pair(avi::lineStatus::coords, c);
      return std::make_pair(avi::lineStatus::garbage, c);
    }
    try {
      auto curVal = std::stod(std::string{first, second});
      if (i > 2) {
        c[0] = c[1]; c[1] = c[2]; c[2] = curVal; 
      } else {
        c[i] = curVal;
      }
    } catch (const std::invalid_argument &) {
      return std::make_pair(avi::lineStatus::garbage, c);
    } catch (const std::out_of_range &) {
      return std::make_pair(avi::lineStatus::garbage, c);
    }
  }
  return std::make_pair(avi::lineStatus::coords, c);
}

std::pair<avi::lineStatus, std::array<avi::Coords, 2>>
avi::getCoordDisplaced(const std::string &line) {
  std::array<avi::Coords, 2> c;
  auto first = std::find_if(begin(line), end(line),
                            [](int ch) { return !std::isspace(ch); });
  if (first == std::end(line))
    return std::make_pair(avi::lineStatus::garbage,
                          c); // possibly blank line
  auto second =
      std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  std::string word{first, second};
  if (word == "ITEM:") {
    first = std::find_if(second, end(line), [](int ch) { return !std::isspace(ch); });
    second = std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    word = std::string{first, second};
    if (word == "ENTRIES") {
      return std::make_pair(avi::lineStatus::frameBorder, c);
    } else {
      return std::make_pair(avi::lineStatus::garbage, c);
    }
  }
  for (auto j = 0; j < 2; ++j)
    for (auto i = 0; i < 3; ++i) {
      first = std::find_if(second, end(line),
                           [](int ch) { return !std::isspace(ch); });
      second = std::find_if(first, end(line),
                            [](int ch) { return std::isspace(ch); });
      if (first >= second)
        return std::make_pair(avi::lineStatus::garbage, c);
      try {
        c[j][i] = std::stod(std::string{first, second});
      } catch (const std::invalid_argument &) {
        return std::make_pair(avi::lineStatus::garbage, c);
      } catch (const std::out_of_range &) {
        return std::make_pair(avi::lineStatus::garbage, c);
      }
    }
  return std::make_pair(avi::lineStatus::coords, c);
}
