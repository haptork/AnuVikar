#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>

#include <helper.hpp>
#include <infoReader.hpp>
#include <printJson.hpp>
#include <results.hpp>
#include <xyz2defects.hpp>

std::tuple<avi::ErrorStatus, avi::InputInfo, avi::ExtraInfo> avi::cookInfos(std::string xyzfileName, avi::infoFrom infoStatus) {
  std::string infileName, tag;
  auto es = avi::ErrorStatus::noError;
  auto info = InputInfo{};
  auto extraInfo = ExtraInfo{};
  auto isInfo = true;
  if (infoStatus == infoFrom::stdinOrFile || infoStatus == infoFrom::file) {
    std::tie(infileName, tag) = avi::getInfileFromXyzfile(xyzfileName);
    if (!infileName.empty()) {
      bool status;
      avi::XyzFileType sc {avi::XyzFileType::generic};
      std::tie(sc, status) = avi::getSimulationCode(infileName);
      if (status) {
        std::tie(info, extraInfo, isInfo) =
        (sc == avi::XyzFileType::parcasWithStdHeader)
            ? avi::extractInfoParcas(infileName, tag)
            : avi::extractInfoLammps(infileName, tag);
        if (isInfo) {
          info.xyzFilePath = xyzfileName;
          info.xyzFileType = sc;
          return std::make_tuple(avi::ErrorStatus::noError, info, extraInfo);
        } else {
          es = avi::ErrorStatus::InputFileincomplete;
        }
      } else {
        es = avi::ErrorStatus::unknownSimulator;
      }
    } else {
      es = avi::ErrorStatus::inputFileMissing;
    }
  }
  if (infoStatus == infoFrom::stdinOrFile || infoStatus == infoFrom::stdin) {
    std::tie(info, extraInfo, isInfo) = avi::infoFromStdIn();
    if (!isInfo) {
      if (es == avi::ErrorStatus::noError) {
        es = avi::ErrorStatus::unknownError;
      } else {
        es = avi::ErrorStatus::noError;
      } 
    }
    return std::make_tuple(avi::ErrorStatus::noError, info, extraInfo);
  }
  return std::make_tuple(es, info, extraInfo);
  // if (isDefaultInfo) Logger::inst().log_info("Found input file " + infileName);
}

std::array<std::string, 2> separateDirAndFile(std::string path) {
  std::size_t dirPos = path.find_last_of("/");
  std::string dir{""};
  std::string file = path;
  if (dirPos != std::string::npos) {
    dir = path.substr(0, dirPos);
    file = path.substr(dirPos + 1);
  }
  return std::array<std::string, 2>{{dir, file}};
}

std::array<std::string, 2> separateFileAndExt(std::string path) {
  std::size_t extPos = path.find_last_of(".");
  std::string fname;
  if (extPos != std::string::npos)
    fname = path.substr(0, extPos);
  else
    return std::array<std::string, 2>{{path, fname}};
  std::string ext = path.substr(extPos + 1);
  return std::array<std::string, 2>{{fname, ext}};
}

std::pair<std::string, std::string>
avi::getInfileFromXyzfile(std::string xyzfile) {
  auto dirFile = separateDirAndFile(xyzfile);
  auto fnameExt = separateFileAndExt(dirFile[1]);
  auto addum = (dirFile[0].length() > 0) ? "/" : "";
  auto origFnameExt = dirFile[1];
  auto origFname = fnameExt[0];
  avi::replaceStr(dirFile[1], "fpos", "md");
  avi::replaceStr(fnameExt[0], "fpos", "md");
  auto filePaths = std::array<std::string, 6>{
      {dirFile[0] + addum + fnameExt[0] + ".in",
       dirFile[0] + addum + dirFile[1] + ".in",
       dirFile[0] + addum + origFnameExt + ".in",
       dirFile[0] + addum + origFname + ".in",
       dirFile[0] + addum + "common_input.in", "common_input.in"}};
  for (const auto &filePath : filePaths) {
    std::ifstream infile(filePath);
    if (!infile.bad() && infile.is_open()) {
      infile.close();
      return std::make_pair(filePath, origFnameExt);
    }
    if (infile.is_open()) infile.close();
  }
  return std::make_pair(std::string{}, origFnameExt);
}

// extract information from input file
std::tuple<avi::InputInfo, avi::ExtraInfo, bool>
avi::extractInfoLammps(std::string fpath, std::string ftag) {
  std::ifstream infile(fpath);
  avi::InputInfo mainInfo;
  avi::ExtraInfo extraInfo;
  if (infile.bad() || !infile.is_open()) {
    return std::make_tuple(mainInfo, extraInfo, false);
  }
  std::string line;
  int xyzrec{0};
  bool careAboutAngles = false;
  auto isOrigin = 0;
  avi::Coords velocity{{0.0, 0.0, 0.0}};
  while (std::getline(infile, line)) {
    auto eq =
        std::find_if(begin(line), end(line), [](int ch) { return ch == '='; });
    if (eq != std::end(line)) {
      auto cmd = trim(std::string{std::begin(line), eq});
      ++eq;
      auto val = removeParcasComments(trim(std::string{eq, std::end(line)}));
      if (cmd == "substrate") {
        extraInfo.substrate = val;
      } else if (cmd == "isPerfect") {
        mainInfo.isPerfect = std::stoi(val);
      } else if (cmd == "boxSizeX") {
        mainInfo.boxSizeX = std::stod(val);
      } else if (cmd == "boxSizeY") {
        mainInfo.boxSizeY = std::stod(val);
      } else if (cmd == "boxSizeZ") {
        mainInfo.boxSizeZ = std::stod(val);
      } else if (cmd == "boxSize") {
        mainInfo.boxSizeX = std::stod(val);
        mainInfo.boxSizeY = std::stod(val);
        mainInfo.boxSizeZ = std::stod(val);
      } else if (cmd == "latticeConstant") {
        mainInfo.latticeConst = std::stod(val);
      } else if (cmd == "ncellX") {
        mainInfo.ncellX = std::stoi(val);
      } else if (cmd == "ncellY") {
        mainInfo.ncellY = std::stoi(val);
      } else if (cmd == "ncellZ") {
        mainInfo.ncellZ = std::stoi(val);
      } else if (cmd == "ncell") {
        mainInfo.ncellX = std::stoi(val);
        if(mainInfo.ncellY < 0) mainInfo.ncellY = mainInfo.ncellX;
        if(mainInfo.ncellZ < 0) mainInfo.ncellZ = mainInfo.ncellX;
      } else if (cmd == "recen") {
        extraInfo.energy = std::stod(val);
      } else if (cmd == "xrec") {
        extraInfo.xrec = std::stod(val);
        xyzrec++;
      } else if (cmd == "yrec") {
        extraInfo.yrec = std::stod(val);
        xyzrec++;
      } else if (cmd == "zrec") {
        extraInfo.zrec = std::stod(val);
        xyzrec++;
      } else if (cmd == "vx") {
        velocity[0] = std::stod(val);
        careAboutAngles = true;
      } else if (cmd == "vy") {
        velocity[1] = std::stod(val);
        careAboutAngles = true;
      } else if (cmd == "vz") {
        velocity[2] = std::stod(val);
        careAboutAngles = true;
      } else if (cmd == "offset") {
        mainInfo.originX = std::stod(val);
        mainInfo.originY = std::stod(val);
        mainInfo.originZ = std::stod(val);
        isOrigin = 3;
      } else if (cmd == "offsetX") {
        mainInfo.originX = std::stod(val);
        isOrigin++;
      } else if (cmd == "offsetY") {
        mainInfo.originY = std::stod(val);
        isOrigin++;
      } else if (cmd == "offsetZ") {
        mainInfo.originZ = std::stod(val);
        isOrigin++;
      } else if (cmd == "temp") {
        mainInfo.temperature = std::stod(val);
      } else if (cmd == "es") {
        extraInfo.es = (val == "true" || val == "yes");
      } else if (cmd == "offsetToUse") {
        if (val == "0") mainInfo.originType = 0;
        if (val == "1") mainInfo.originType = 1;
        if (val == "2") mainInfo.originType = 2;
      } else if (cmd == "author") {
        extraInfo.author = val;
      } else if (cmd == "potentialUsed") {
        extraInfo.potentialUsed = val;
      } else if (cmd == "xColumn") {
        mainInfo.xyzColumnStart = std::stoi(val);
      } else if (cmd == "structure") {
        mainInfo.structure = val;
      } else if (cmd == "frameStart") {
        mainInfo.frameStart = std::stoi(val);
      } else if (cmd == "frameEnd") {
        mainInfo.frameEnd = std::stoi(val);
      } else if (cmd == "framePeriod") {
        mainInfo.framePeriod = std::stoi(val);
      } else if (cmd == "extraColumn") {
        if (mainInfo.extraColumnStart == -1) { // start only
          mainInfo.extraColumnStart = std::stoi(val);
          mainInfo.extraColumnEnd = mainInfo.extraColumnStart;
        } else {
          mainInfo.extraColumnEnd = std::stoi(val);;
        }
      }
      /*
            } else if (cmd == "latConstToUse") {
              if (val == "0") mainInfo.latConstType = 0;
              if (val == "1") mainInfo.latConstType = 1;
              if (val == "2") mainInfo.latConstType = 2;
      */
    }
  }
  infile.close();
  if (isOrigin == 3) mainInfo.originType = 0;
  extraInfo.isPkaGiven = (xyzrec >= 3) ? true : false;
  if (mainInfo.latticeConst < 0.0) {
    return std::make_tuple(mainInfo, extraInfo, false);
  }
  if (mainInfo.ncellX > 0 && mainInfo.isPerfect) {
    mainInfo.boxSizeX = mainInfo.latticeConst * mainInfo.ncellX;
    mainInfo.boxSizeY = mainInfo.latticeConst * mainInfo.ncellY;
    mainInfo.boxSizeZ = mainInfo.latticeConst * mainInfo.ncellZ;
  }
  if (mainInfo.xyzColumnStart == -1 && mainInfo.extraColumnStart > -1) {
    return std::make_tuple(mainInfo, extraInfo, false);
  }
  extraInfo.infile = fpath;
  extraInfo.rectheta =
      (careAboutAngles) ? std::atan(velocity[1] / velocity[0]) : 0.0;
  extraInfo.recphi =
      (careAboutAngles) ? std::atan(velocity[2] / velocity[0]) : 0.0;
  return std::make_tuple(mainInfo, extraInfo, true);
}

std::string readStr(std::string msg) {
  std::string buffer;
  std::cout << msg;
  std::getline(std::cin, buffer);
  return buffer;
}

template <size_t N>
std::pair<bool, std::array<double,  N>> readAr(std::string msg) {
  auto line = readStr(msg);
  std::array<double, N> res;
  if (line.empty()) return std::make_pair(false, res);
  auto first = std::find_if(begin(line), end(line),
                            [](int ch) { return !std::isspace(ch); });
  if (first == std::end(line)) return std::make_pair(false, res); // possibly blank line
  auto second =
      std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
  try {
    res[0] = std::stod(std::string{first, second});
    for (int i = 1; i < N; ++i) res[i] = res[0];
  } catch (const std::invalid_argument &) {
    return std::make_pair(false, res);
  } catch (const std::out_of_range &) {
    return std::make_pair(false, res);
  }
  for (int i = 1; i < N; ++i) {
    first = std::find_if(second, end(line),
                         [](int ch) { return !std::isspace(ch); });
    second =
        std::find_if(first, end(line), [](int ch) { return std::isspace(ch); });
    if (first >= second) return std::make_pair(true, res);
    try {
      res[i] = std::stod(std::string{first, second});
    } catch (const std::invalid_argument &) {
      return std::make_pair(false, res);
    } catch (const std::out_of_range &) {
      return std::make_pair(false, res);
    }
  }
  return std::make_pair(true, res);
}

double readDouble(std::string msg) {
  auto buffer = readStr(msg);
  if (buffer.empty()) return 0.0;
  double res = 0.0;
  try {
      res = std::stod(buffer);
  } catch (const std::invalid_argument &) {
    std::cout << "Input is not valid number.";
    buffer.clear();
  } catch (const std::out_of_range &) {
    std::cout << "Input is not valid number.";
    buffer.clear();
  }
  return res;
}

int readInt(std::string msg) {
  auto buffer = readStr(msg);
  if (buffer.empty()) return 0;
  int res = 0;
  try {
      res = std::stoi(buffer);
  } catch (const std::invalid_argument &) {
    std::cout << "Input is not valid number.";
    buffer.clear();
  } catch (const std::out_of_range &) {
    std::cout << "Input is not valid number.";
    buffer.clear();
  }
  return res;

}

std::tuple<avi::InputInfo, avi::ExtraInfo, bool> avi::infoFromStdIn() {
  avi::InputInfo mainInfo;
  avi::ExtraInfo extraInfo;
  std::string buffer;
  std::cout << "\rInput file missing, please provide inputs here: \n"; 
  while (buffer.size() == 0 || mainInfo.latticeConst < 0.0) {
    std::cout << "Lattice constant: ";
    std::getline(std::cin, buffer);
    try {
      mainInfo.latticeConst = std::stod(buffer);
    } catch (const std::invalid_argument &) {
      std::cout << "Input is not a valid number.";
      buffer.clear();
    } catch (const std::out_of_range &) {
      std::cout << "Input is not a valid number.";
      buffer.clear();
    }
  }
  std::cout << "Optional parameters (Press return / enter to continue with default): \n";
  auto origin = readAr<3>("offset used for simulation: ");
  mainInfo.originType = origin.first ? 0 : 1;
  if (origin.first) {
    mainInfo.originX = origin.second[0];
    mainInfo.originY = origin.second[1];
    mainInfo.originZ = origin.second[2];
  }
  extraInfo.substrate = readStr("substrate symbol (e.g.: W/Fe): ");
  mainInfo.structure = readStr("structure (bcc / fcc (default: bcc)): ");
  extraInfo.energy = readDouble("PKA energy (in keV): ");
  mainInfo.temperature = readDouble("Temperature (K): ");
  extraInfo.author = readStr("Author name: ");
  extraInfo.potentialUsed = readStr("Potential used: ");
  extraInfo.es = readStr("Is electronic stopping used (y/n): ") == "y" ? true : false;
  auto xyzFormat = readInt("xyz file formate (1-generic-xyz (default), 2-lammps, 3-parcas, 4-cascadesDb 5-displaced): ");
  std::vector<avi::XyzFileType> codes{
    avi::XyzFileType::generic,
    avi::XyzFileType::lammpsWithStdHeader,
    avi::XyzFileType::parcasWithStdHeader,
    avi::XyzFileType::cascadesDbLikeCols,
    avi::XyzFileType::lammpsDisplacedCompute
  };
  mainInfo.xyzFileType = (xyzFormat < codes.size() && xyzFormat > 0) ? codes[xyzFormat - 1] : codes[0];
  auto xyzCol = readInt("column number for x-coordinate (subsequent columns will be taken as y & z) (default: auto): ");
  mainInfo.xyzColumnStart = (xyzFormat > 0) ? xyzCol : -1;
  auto ecInfo = readAr<2>("extra columns (default: none): ");
  if (ecInfo.first) {
    mainInfo.extraColumnStart = ecInfo.second[0];
    mainInfo.extraColumnEnd = ecInfo.second[1];
  }
  return std::make_tuple(mainInfo, extraInfo, true);
}

// extract information from input file
std::tuple<avi::InputInfo, avi::ExtraInfo, bool>
avi::extractInfoParcas(std::string fpath, std::string ftag) {
  std::ifstream infile(fpath);
  avi::InputInfo mainInfo;
  avi::ExtraInfo extraInfo;
  if (!infile.bad() && infile.is_open()) {
    std::string line;
    auto count = 0;
    while (std::getline(infile, line)) {
      auto i = line.find('=');
      if (i == 9) {
        auto cmd = trim(line.substr(0, 9));
        if (cmd == "substrate") {
          extraInfo.substrate = removeParcasComments(line.substr(10));
          count++;
        } else if (cmd == "box(1)") {
          mainInfo.boxSizeX = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "box(2)") {
          mainInfo.boxSizeY = std::stod(removeParcasComments(line.substr(10)));
          count++;
         } else if (cmd == "box(3)") {
          mainInfo.boxSizeZ = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "ncell(1)") {
          mainInfo.ncellX = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "ncell(2)") {
          mainInfo.ncellY = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "ncell(3)") {
          mainInfo.ncellZ = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "recen") {
          extraInfo.energy =
              (std::stod(removeParcasComments(line.substr(10)))) / 1000.0;
          count++;
        } else if (cmd == "xrec") {
          extraInfo.xrec = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "yrec") {
          extraInfo.yrec = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "zrec") {
          extraInfo.zrec = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "rectheta") {
          extraInfo.rectheta = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "recphi") {
          extraInfo.recphi = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "offset(1)") {
          mainInfo.originX = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "offset(2)") {
          mainInfo.originY = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "offset(3)") {
          mainInfo.originZ = std::stod(removeParcasComments(line.substr(10)));
          count++;
        } else if (cmd == "temp") {
          mainInfo.temperature =
              std::stod(removeParcasComments(line.substr(10)));
          count++;
        }
      }
    }
    mainInfo.structure = "bcc"; // TODO: extend for fcc, without assumption
    infile.close();
    if (count == 17) { // got all the info
      mainInfo.latticeConst = mainInfo.boxSizeX / mainInfo.ncellX;
      extraInfo.isPkaGiven = true;
      extraInfo.infile = fpath;
      return std::make_tuple(mainInfo, extraInfo, true);
    }
    return std::make_tuple(mainInfo, extraInfo, false);
  }
  return std::make_tuple(mainInfo, extraInfo, false);
}

std::pair<avi::XyzFileType, bool>
avi::getSimulationCode(std::string fname) {
  std::ifstream infile(fname);
  if (!infile.bad() && infile.is_open()) {
    std::string line;
    auto lineNo = 0;
    constexpr auto maxLinesToLook = 10;
    std::vector<std::string> keyWords{ "CASCADESDBLIKECOLS", "PARCAS",
                                      "LAMMPS-XYZ", "LAMMPS-DISP", "XYZ",};
    std::vector<avi::XyzFileType> codes{
        avi::XyzFileType::cascadesDbLikeCols,
        avi::XyzFileType::parcasWithStdHeader,
        avi::XyzFileType::lammpsWithStdHeader,
        avi::XyzFileType::lammpsDisplacedCompute,
        avi::XyzFileType::generic,
        };
    while (std::getline(infile, line) && lineNo++ < maxLinesToLook) {
      for (size_t i = 0; i < keyWords.size(); i++) {
        auto pos = line.find(keyWords[i]);
        if (pos != std::string::npos) return std::make_pair(codes[i], true);
      }
    }
  }
  return std::make_pair(avi::XyzFileType::generic, false);
}
