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

/*
* Get infos from input file if input file is given else uses defaultInfo if given.
* 
*/
std::pair<avi::ErrorStatus,int> avi::processFileTimeCmd(std::string xyzfileName,
                                            std::ostream &outfile,
                                            const avi::Config &config, int id, avi::InputInfo info, avi::ExtraInfo extraInfo, avi::infoFrom infoStatus) {
  auto es = avi::ErrorStatus::noError;
  auto isInfo = true;
  if (infoStatus != avi::infoFrom::cooked) {
    std::tie(es, info, extraInfo) = cookInfos(xyzfileName, infoStatus);
    if (es != avi::ErrorStatus::noError) return std::make_pair(es, 0);
  }
  avi::frameStatus fs = avi::frameStatus::prelude;
  std::ifstream xyzfile{info.xyzFilePath};
  if (xyzfile.bad() || !xyzfile.is_open()) return std::make_pair(avi::ErrorStatus::xyzFileOpenError, 0);
  auto success = 0;
  auto frameCount = 0;
  while (true) {
    extraInfo.simulationTime = frameCount;
    extraInfo.id = std::to_string(id + success + 1);
    if (frameCount < info.frameStart || frameCount % info.framePeriod != 0) {
      auto res = avi::skipFrame(info, extraInfo, xyzfile, fs);
      frameCount++;
      if (info.frameEnd > 0 && frameCount >= info.frameEnd || res == avi::xyzFileStatus::eof) break;
      continue;
    }
    auto res = avi::processTimeFile(info, extraInfo, config, xyzfile, fs, outfile, success == 0);
    if (res.second != avi::ErrorStatus::noError) {
      if (config.allFrames) std::cerr << "\nError: " << errToStr(res.second) << " in frame " << frameCount << " of file " << xyzfileName << " input-file: " << extraInfo.infile << '\n' << std::flush;
      else std::cerr << "\nError: " << errToStr(res.second) << " : file- " << xyzfileName << " input-file: " << extraInfo.infile << '\n' << std::flush;
      Logger::inst().log_info("Error processing frame no. " + std::to_string(frameCount) +" in file \"" + info.xyzFilePath + "\"");
    } else {
      ++success;
      if (config.allFrames) {
        if (success >= 2) std::cout << "\r" << frameCount << " step processed successfully." << std::flush;
        Logger::inst().log_info("Finished processing frame no. " + std::to_string(success) +" in file \"" + xyzfileName + "\"");
      }
    }
    frameCount++;
    if ((info.frameEnd > 0 && frameCount >= info.frameEnd) || res.first == avi::xyzFileStatus::eof) break;
  }
  xyzfile.close();
  if (success > 0) return std::make_pair(avi::ErrorStatus::noError, success);
  return std::make_pair(avi::ErrorStatus::unknownError, 0);
}

inline int countDefects(const avi::DefectVecT &defects) {
  return std::count_if(begin(defects), end(defects), [](const auto &x) {
    return avi::DefectTWrap::isSurviving(x) && avi::DefectTWrap::isInterstitial(x);
  });
}

std::pair<avi::xyzFileStatus, avi::ErrorStatus> 
                          avi::processTimeFile(avi::InputInfo &info,
                                     avi::ExtraInfo &extraInfo,
                                     const avi::Config &config, std::istream &infile, avi::frameStatus &fs, std::ostream &outfile, bool isFirst) {
  auto res = avi::resultsT{};
  //res.err = avi::ErrorStatus::noError;
  avi::xyzFileStatus fl;
  std::cout << "her\n";
  std::tie(fl, res.err, res.defects, res.coDefects, res.extraCols) = 
      (info.xyzFileType == avi::XyzFileType::lammpsDisplacedCompute)
          ? avi::displaced2defectsTime(info, extraInfo, config, infile, fs)
          : avi::xyz2defectsTime(info, extraInfo, config, infile, fs);
  if (res.err != avi::ErrorStatus::noError) return std::make_pair(fl, res.err);
  if (config.onlyDefects) {
    if (!isFirst) outfile << "\n,";
    res.nDefects = countDefects(res.defects);
    res.nSia = res.nDefects;
    res.nVac = res.nDefects;
    avi::printJson(outfile, info, extraInfo, res);
    return std::make_pair(fl, res.err);
  }
  avi::Coords box{{info.boxSizeX, info.boxSizeY, info.boxSizeZ}};
  res.defects = avi::groupDefects(std::move(res.defects), info.latticeConst, box);
  //std::cout << "\nec size at reader: " << res.extraCols.size() << '\n';
  auto clusterSizeMap = avi::clusterSizes(res.defects);
  filterZeroClusters(res.defects, clusterSizeMap,
                     config.filterZeroSizeClusters);
  avi::ignoreSmallClusters(res.defects, clusterSizeMap);
  res.clusters = avi::clusterMapping(res.defects);
  res.clustersIV = avi::clusterIVType(res.clusters, clusterSizeMap);
  if (config.isFindClusterFeatures)
    res.feats = avi::clusterFeatures(res.defects, res.clusters,
                                          clusterSizeMap, info.latticeConst, box);
  std::tie(res.nSia, res.nVac, res.inClusterFractionI, res.inClusterFractionV) =
      avi::getNDefectsAndClusterFractions(res.defects);
  res.nDefects = std::min(res.nSia, res.nVac);
  std::tie(res.maxClusterSizeI, res.maxClusterSizeV) =
      avi::getMaxClusterSizes(clusterSizeMap, res.clusters);
  res.nClusters = res.clusters.size();
  if (res.err == avi::ErrorStatus::noError) {
    if (!isFirst) outfile << "\n,";
    avi::printJson(outfile, info, extraInfo, res);
  }
  return std::make_pair(fl, res.err);
}
