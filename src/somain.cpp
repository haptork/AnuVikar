#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <helper.hpp>
#include <logger.hpp>
#include <printJson.hpp>
#include <reader.hpp>
#include <results.hpp>
#include <infoReader.hpp>
#include <xyz2defects.hpp>

struct InfoPyInput {
  int ncellX;
  int ncellY;
  int ncellZ;
  double boxSizeX;
  double boxSizeY;
  double boxSizeZ;
  double latticeConst;
  double originX;
  double originY;
  double originZ;
  int originType;
  double temperature;
  int xyzColumnStart;
  int extraColumnStart;
  int extraColumnEnd;
  int frameStart{-1};
  int frameEnd{-1};
  int framePeriod{1};
  bool isPerfect;
  const char *xyzFileType;
  const char *xyzFilePath;
  const char *structure;
};
struct InfoPyExtraInput {
  double energy;
  double simulationTime;
  // for distribution around PKA
  bool isPkaGiven;
  bool es;
  double xrec{0.0};
  double yrec{0.0};
  double zrec{0.0};
  double rectheta{0.0};
  double recphi{0.0};
  // extra
  const char *id;
  const char *substrate;
  const char *infile;
  const char *tags;
  const char *potentialUsed;
  const char *author;
};

struct PyConfig {
  bool onlyDefects{false};
  bool isFindDistribAroundPKA{true};
  bool isFindClusterFeatures{true};
  bool filterZeroSizeClusters{false};
  bool isIgnoreBoundaryDefects{false};
  bool isAddThresholdDefects{true};
  bool safeRunChecks{true};
  int allFrames{0};
  double thresholdFactor{0.345};
  double extraDefectsSafetyFactor{50.0};
  int logMode{avi::LogMode::warning | avi::LogMode::error};
  const char *logFilePath;
  const char *outputJSONFilePath;
};

auto pyConfigToCppConfig(const PyConfig &pyConfig) {
  auto config = avi::Config{};
  config.allFrames = pyConfig.allFrames;
  config.onlyDefects = pyConfig.onlyDefects;
  config.isFindClusterFeatures = pyConfig.isFindClusterFeatures;
  config.isFindDistribAroundPKA = pyConfig.isFindDistribAroundPKA;
  config.filterZeroSizeClusters = pyConfig.filterZeroSizeClusters;
  config.isIgnoreBoundaryDefects = pyConfig.isIgnoreBoundaryDefects;
  config.isAddThresholdDefects = pyConfig.isAddThresholdDefects;
  config.safeRunChecks = pyConfig.safeRunChecks;
  config.thresholdFactor = pyConfig.thresholdFactor;
  config.extraDefectsSafetyFactor = pyConfig.extraDefectsSafetyFactor;
  config.logMode = pyConfig.logMode;
  config.logFilePath = std::string{pyConfig.logFilePath};
  config.outputJSONFilePath = std::string{pyConfig.outputJSONFilePath};
  return config;
}

auto pyInfoToCppInfo(const InfoPyInput &pyinput,
                     const InfoPyExtraInput &pyextra) {
  avi::InputInfo input;
  input.ncellX = pyinput.ncellX;
  input.ncellY = pyinput.ncellY;
  input.ncellZ = pyinput.ncellZ;
  input.boxSizeX = pyinput.boxSizeX;
  input.boxSizeY = pyinput.boxSizeY;
  input.boxSizeZ = pyinput.boxSizeZ;
  input.latticeConst = pyinput.latticeConst;
  input.originX = pyinput.originX;
  input.originY = pyinput.originY;
  input.originZ = pyinput.originZ;
  input.originType = pyinput.originType;
  input.temperature = pyinput.temperature;
  input.xyzColumnStart = pyinput.xyzColumnStart;
  input.extraColumnStart = pyinput.extraColumnStart;
  if (pyinput.extraColumnStart != -1 && pyinput.extraColumnEnd == -1) {
      input.extraColumnEnd = pyinput.extraColumnStart;
  } else {
      input.extraColumnEnd = pyinput.extraColumnEnd;
  }
  input.frameStart = pyinput.frameStart;
  input.frameEnd = pyinput.frameEnd;
  input.framePeriod = pyinput.framePeriod;
  input.isPerfect = pyinput.isPerfect;
  input.xyzFilePath = std::string{pyinput.xyzFilePath};
  input.structure = std::string{pyinput.structure};
  std::string simCodeStr = std::string{pyinput.xyzFileType};
  std::transform(simCodeStr.begin(), simCodeStr.end(), simCodeStr.begin(),
                 [](unsigned char c) { return std::toupper(c); } // correct
  );
  std::vector<std::string> keyWords{"GENERIC", "CASCADESDBLIKECOLS", "PARCAS",
                                    "LAMMPS-XYZ", "LAMMPS-DISP"};
  std::vector<avi::XyzFileType> codes{
      avi::XyzFileType::generic,
      avi::XyzFileType::cascadesDbLikeCols,
      avi::XyzFileType::parcasWithStdHeader,
      avi::XyzFileType::lammpsWithStdHeader,
      avi::XyzFileType::lammpsDisplacedCompute};
  for (size_t i = 0; i < keyWords.size(); i++) {
    if (simCodeStr == keyWords[i]) {
      input.xyzFileType = codes[i];
      break;
    }
  }
  avi::ExtraInfo extra;
  extra.energy = pyextra.energy;
  extra.simulationTime = pyextra.simulationTime;
  extra.isPkaGiven = pyextra.isPkaGiven;
  extra.es = pyextra.es;
  extra.xrec = pyextra.xrec;
  extra.yrec = pyextra.yrec;
  extra.zrec = pyextra.zrec;
  extra.rectheta = pyextra.rectheta;
  extra.recphi = pyextra.recphi;
  extra.id = std::string{pyextra.id};
  extra.substrate = std::string{pyextra.substrate};
  extra.infile = std::string{pyextra.infile};
  extra.tags = std::string{pyextra.tags};
  extra.potentialUsed = std::string{pyextra.potentialUsed};
  extra.author = std::string{pyextra.author};
  return std::make_tuple(true, input, extra);
}

extern "C" char *pyProcessFile(InfoPyInput pyInfo, InfoPyExtraInput pyExtraInfo,
                               PyConfig pyConfig) {
  using avi::Logger;
  auto config = pyConfigToCppConfig(pyConfig);
  auto info = avi::InputInfo{};
  auto extraInfo = avi::ExtraInfo{};
  auto isSuccess = false;
  std::tie(isSuccess, info, extraInfo) = pyInfoToCppInfo(pyInfo, pyExtraInfo);
  Logger::inst().mode(config.logMode);
  Logger::inst().file(config.logFilePath);
  std::stringstream outfile;
  auto res = avi::resultsT{};
  auto lastErr = avi::ErrorStatus::noError;
  auto success = 0;
  if (config.allFrames) outfile << "[";
  if (!isSuccess) {
    res.err = avi::ErrorStatus::unknownSimulator;
  } else {
    Logger::inst().log_info("Started Processing file \"" + info.xyzFilePath +
                            " (" + extraInfo.infile + ") " + "\"");
    avi::frameStatus fs = avi::frameStatus::prelude;
    std::ifstream xyzfile{info.xyzFilePath};
    if (xyzfile.bad() || !xyzfile.is_open()) {
      res.err = avi::ErrorStatus::xyzFileOpenError;
    }
    auto initId = extraInfo.id;
    auto origSimTime = extraInfo.simulationTime;
    auto frameCount = 0;
    while (true) {
      extraInfo.simulationTime = frameCount + origSimTime;
      extraInfo.id = initId + "_" + std::to_string(success + 1);
      if (frameCount < info.frameStart || frameCount % info.framePeriod != 0) {
        auto res = avi::skipFrame(info, extraInfo, xyzfile, fs);
        frameCount++;
        if (info.frameEnd > 0 && frameCount >= info.frameEnd || res == avi::xyzFileStatus::eof) break;
        continue;
      }
      auto res = avi::processTimeFile(info, extraInfo, config, xyzfile, fs, outfile, success == 0);
      if (res.second != avi::ErrorStatus::noError) {
        Logger::inst().log_info("Error processing frame no." + std::to_string(frameCount) +" in file \"" + info.xyzFilePath + "\"");
        lastErr = res.second;
      } else {
        ++success;
        if (config.allFrames) 
        Logger::inst().log_info("Finished processing frame no. " + std::to_string(frameCount) +" in file \"" + info.xyzFilePath + "\"");
      }
      frameCount++;
      if ((info.frameEnd > 0 && frameCount >= info.frameEnd) || res.first == avi::xyzFileStatus::eof) break;
    }
    xyzfile.close();
    Logger::inst().log_info("Finished Processing");
  }
  if (success == 0) {
    res.err = lastErr;
    avi::printJson(outfile, info, extraInfo, res);
  }
  if (config.allFrames) outfile << "]";
  std::string str = outfile.str();
  char *writable = (char *)malloc(sizeof(char) * (str.size() + 1));
  std::copy(str.begin(), str.end(), writable);
  writable[str.size()] = '\0';
  return writable;
}

extern "C" char *pyProcessFileWoInfo(char *xyzfile, PyConfig pyConfig) {
  using avi::Logger;
  auto config = pyConfigToCppConfig(pyConfig);
  auto xyzfileStr = std::string(xyzfile);
  Logger::inst().mode(config.logMode);
  Logger::inst().file(config.logFilePath);
  Logger::inst().log_info("Started Processing file \"" + xyzfileStr + "\"");
  std::stringstream outfile;
  avi::InputInfo info;
  avi::ExtraInfo extraInfo;
  avi::infoFrom infoStatus = avi::infoFrom::file;
  auto res = avi::processFileTimeCmd(xyzfileStr, outfile, config, 0, info, extraInfo, infoStatus);
  Logger::inst().log_info("Finished Processing");
  std::string str = outfile.str();
  char *writable = (char *)malloc(
      sizeof(char) * (str.size() + 1)); // new char[str.size() + 1];
  std::copy(str.begin(), str.end(), writable);
  writable[str.size()] = '\0';
  return writable;
}

extern "C" void dalloc(void *x) { free(x); }
