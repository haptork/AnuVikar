#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>

#include <helper.hpp>
#include <infoReader.hpp>
#include <logger.hpp>
#include <printJson.hpp>
#include <reader.hpp>
#include <results.hpp>
#include <xyz2defects.hpp>

auto filterZeroClusters(avi::DefectVecT &defects,
                        avi::ClusterSizeMapT &clusterSize, bool isFilter) {
  using namespace avi::DefectTWrap;
  if (avi::Logger::inst().mode() & avi::LogMode::debug) {
    for (auto it : clusterSize) {
      if (it.second.surviving == 0) {
        avi::Logger::inst().log_debug(
            "Found cluster with zero size. (id, total defects): " +
            std::to_string(it.first) + ", " + std::to_string(it.second.all));
      }
    }
  }
  if (!isFilter) return;
  auto initSize = defects.size();
  defects.erase(
      std::remove_if(defects.begin(), defects.end(),
                     [&clusterSize](const auto &defect) {
                       return clusterSize[clusterId(defect)].surviving == 0;
                     }),
      defects.end());
}

std::pair<avi::ErrorStatus,int> avi::processFileTimeCmd(std::string xyzfileName,
                                            std::ostream &outfile,
                                            const avi::Config &config, int id, const avi::InputInfo &defaultInfo, const avi::ExtraInfo &defaultExtraInfo, bool isDefaultInfo) {
  std::string infileName, tag;
  std::tie(infileName, tag) = avi::getInfileFromXyzfile(xyzfileName);
  //if (infileName.empty()) return std::make_pair(avi::ErrorStatus::inputFileMissing, 0);
  avi::XyzFileType sc {avi::XyzFileType::generic};
  avi::InputInfo info;
  avi::ExtraInfo extraInfo;
  bool isInfo;
  if (infileName.empty()) {
    if (!isDefaultInfo) return std::make_pair(avi::ErrorStatus::inputFileMissing, 0);
    info = defaultInfo;
    extraInfo = defaultExtraInfo;
    isInfo = isDefaultInfo;
    sc = info.xyzFileType;
  } else {
    bool status;
    std::tie(sc, status) = avi::getSimulationCode(infileName);
    if (!status) return std::make_pair(avi::ErrorStatus::unknownSimulator, 0);
    std::tie(info, extraInfo, isInfo) =
      (sc == avi::XyzFileType::parcasWithStdHeader)
          ? avi::extractInfoParcas(infileName, tag)
          : avi::extractInfoLammps(infileName, tag);
    if (isDefaultInfo) Logger::inst().log_info("Found input file " + infileName);
  }
  if (!isInfo) return std::make_pair(avi::ErrorStatus::InputFileincomplete, 0);
  info.xyzFileType = sc;
  info.xyzFilePath = xyzfileName;
  avi::frameStatus fs = avi::frameStatus::prelude;
  std::ifstream xyzfile{info.xyzFilePath};
  if (xyzfile.bad() || !xyzfile.is_open()) return std::make_pair(avi::ErrorStatus::xyzFileReadError, 0);
  auto success = 0;
  auto frameCount = 0;
  while (true) {
    extraInfo.simulationTime = success + 1;
    extraInfo.id = std::to_string(id + success + 1);
    auto res = avi::processTimeFile(info, extraInfo, config, xyzfile, fs, outfile, success == 0);
    frameCount++;
    if (res.second != avi::ErrorStatus::noError) {
      if (config.allFrames) std::cerr << "\nError: " << errToStr(res.second) << " in frame " << frameCount << " of file " << xyzfileName << '\n' << std::flush;
      else std::cerr << "\nError: " << errToStr(res.second) << " of file " << xyzfileName << '\n' << std::flush;
      Logger::inst().log_info("Error processing" + std::to_string(frameCount) +" frame in file \"" + xyzfileName + "\"");
    } else {
      ++success;
      if (config.allFrames) {
        if (success >= 2) std::cout << "\r" << success << " steps processed successfully." << std::flush;
        Logger::inst().log_info("Finished processing" + std::to_string(success) +" frame in file \"" + xyzfileName + "\"");
      }
    }
    if (res.first == avi::xyzFileStatus::eof) break;
  }
  xyzfile.close();
if (success > 0) return std::make_pair(avi::ErrorStatus::noError, success);
  return std::make_pair(avi::ErrorStatus::unknownError, 0);
}

std::pair<avi::xyzFileStatus, avi::ErrorStatus> 
                          avi::processTimeFile(avi::InputInfo &info,
                                     avi::ExtraInfo &extraInfo,
                                     const avi::Config &config, std::istream &infile, avi::frameStatus &fs, std::ostream &outfile, bool isFirst) {
  auto res = avi::resultsT{};
  //res.err = avi::ErrorStatus::noError;
  avi::xyzFileStatus fl;
  std::tie(fl, res.err, res.defects, res.coDefects, res.extraCols) = 
      (info.xyzFileType == avi::XyzFileType::lammpsDisplacedCompute)
          ? avi::displaced2defectsTime(info, extraInfo, config, infile, fs)
          : avi::xyz2defectsTime(info, extraInfo, config, infile, fs);
  if (res.err != avi::ErrorStatus::noError) return std::make_pair(fl, res.err);
  res.defects = avi::groupDefects(std::move(res.defects), info.latticeConst);
  //std::cout << "\nec size at reader: " << res.extraCols.size() << '\n';
  auto clusterSizeMap = avi::clusterSizes(res.defects);
  filterZeroClusters(res.defects, clusterSizeMap,
                     config.filterZeroSizeClusters);
  avi::ignoreSmallClusters(res.defects, clusterSizeMap);
  res.clusters = avi::clusterMapping(res.defects);
  res.clustersIV = avi::clusterIVType(res.clusters, clusterSizeMap);
  if (config.isFindClusterFeatures)
    res.feats = avi::clusterFeatures(res.defects, res.clusters,
                                          clusterSizeMap, info.latticeConst);
  int nDefects;
  std::tie(res.nDefects, res.inClusterFractionI, res.inClusterFractionV) =
      avi::getNDefectsAndClusterFractions(res.defects);
  std::tie(res.maxClusterSizeI, res.maxClusterSizeV) =
      avi::getMaxClusterSizes(clusterSizeMap, res.clusters);
  res.nClusters = res.clusters.size();
  if (res.err == avi::ErrorStatus::noError) {
    if (!isFirst) outfile << "\n,";
    avi::printJson(outfile, info, extraInfo, res);
  }
  return std::make_pair(fl, res.err);
}
