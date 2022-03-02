#include <cluster2features.hpp>
#include <printJson.hpp>
#include <results.hpp>
#include <sstream>
#include <iomanip>

void printClusterIds(const std::unordered_map<int, std::vector<int>> clusters,
                     std::ostream &outfile) {
  size_t i = 0;
  for (const auto &it : clusters) {
    outfile << '"' << it.first << "\":";
    outfile << "[";
    avi::writeVector(it.second, outfile);
    outfile << "]";
    if (i != clusters.size() - 1) outfile << ", ";
    outfile << "\n";
    ++i;
  }
}

void printClusterIVs(const std::unordered_map<int, int> clusters,
                     std::ostream &outfile) {
  size_t i = 0;
  for (const auto &it : clusters) {
    outfile << '"' << it.first << "\":";
    outfile << it.second;
    if (i != clusters.size() - 1) outfile << ", ";
    outfile << "\n";
    ++i;
  }
}

void printSingleFeat(const avi::featT &feats, std::ostream &outfile) {
  outfile << "\"dist\": ";
  outfile << "[";
  avi::writeStdAr(std::get<0>(feats), outfile);
  outfile << "],\n";
  outfile << "\"angle\": ";
  outfile << "[";
  avi::writeStdAr(std::get<1>(feats), outfile);
  /*
  outfile << "],\n";
  outfile << "\"adjNn2\": ";
  outfile << "[";
  avi::writeStdAr(std::get<2>(feats), outfile);
  */
  outfile << "]\n";
}

void printFeats(const std::unordered_map<int, avi::featT> &feats,
                std::ostream &outfile) {
  size_t count = 0;
  for (const auto &it : feats) {
    outfile << '"' << it.first << "\":";
    outfile << "{";
    printSingleFeat(it.second, outfile);
    outfile << "}";
    if (count++ != feats.size() - 1) outfile << ", ";
    outfile << "\n";
  }
}

auto strSimulationCode(avi::XyzFileType code) {
  return (code == avi::XyzFileType::cascadesDbLikeCols)
             ? "cascadesDbLikeCols"
             : (code == avi::XyzFileType::lammpsWithStdHeader)
                   ? "lammpsWithStdHeader"
                   : (code == avi::XyzFileType::parcasWithStdHeader)
                      ? "parcasWithStdHeader"
                      : (code == avi::XyzFileType::lammpsDisplacedCompute)
                        ? "lammpsDisp"
                        : "generic-XYZ";
}

std::string errorStr(avi::ErrorStatus err) {
  if (err == avi::ErrorStatus::inputFileMissing) {
    return "Could not read input file";
  } else if (err == avi::ErrorStatus::InputFileincomplete) {
    return "Input file doesn't have all the info";
  } else if (err == avi::ErrorStatus::inputFileReadError) {
    return "Input file read error.";
  } else if (err == avi::ErrorStatus::xyzFileMissing) {
    return "Xyz file missing.";
  } else if (err == avi::ErrorStatus::xyzFileReadError) {
    return "Xyz file read err";
  } else if (err == avi::ErrorStatus::unknownSimulator) {
    return "Input file doesn't have LAMMPS/PARCAS/DISPLACED simulation input "
           "type";
  } else if (err == avi::ErrorStatus::xyzFileDefectsProcessingError) {
    return "XYZ file has too many defects or zero atoms";
  } else if (err == avi::ErrorStatus::vacOverflow) {
    return "Too many vacancies";
  } else if (err == avi::ErrorStatus::siaOverflow) {
    return "Too many Sias";
  } else if (err == avi::ErrorStatus::defectOverflow) {
    return "Too many defects";
  } else if (err == avi::ErrorStatus::threshOverflow) {
    return "Too many threshold based defects";
  } else if (err == avi::ErrorStatus::siaVacDiffOverflow) {
    return "Too big difference in sia & vacancy.";
  } else if (err == avi::ErrorStatus::noError) {
    return "";
  }
  return "Unknown error!";
}

void printExtraCols(std::vector<std::vector<double>> ec, std::ostream &outfile) {
  size_t i = 0;
  for (const auto &it : ec) {
    outfile << "[";
    avi::writeVector(it, outfile);
    outfile << "]";
    if (i != ec.size() - 1) outfile << ", ";
    ++i;
  }
}

void avi::resToKeyValue(std::ostream &outfile,
                             const avi::resultsT &res) {
  auto printDefects = [&outfile](const avi::DefectVecT &d) {
    size_t count = 0;
    for (const auto &x : d) {
      outfile << "[" << std::get<0>(x)[0] << ", " << std::get<0>(x)[1] << ", "
              << std::get<0>(x)[2] << ", " << std::get<1>(x) << ", "
              << std::get<2>(x) << ", " << std::get<3>(x) << "]";
      if (count++ < (d.size() - 1)) outfile << ", ";
    }
  };
  outfile << "\"error\":\"" << errorStr(res.err) << "\",\n"
          << "\"n_defects\":" << res.nDefects << ",\n"
          << "\"n_clusters\":" << res.nClusters << ",\n"
          << "\"max_cluster_size_I\":" << res.maxClusterSizeI << ",\n"
          << "\"max_cluster_size_V\":" << res.maxClusterSizeV << ",\n"
          << "\"max_cluster_size\":"
          << std::max(res.maxClusterSizeI, res.maxClusterSizeV) << ",\n"
          << "\"in_cluster_I\":" << res.inClusterFractionI << ",\n"
          << "\"in_cluster_V\":" << res.inClusterFractionV << ",\n"
          << "\"in_cluster\":"
          << (res.inClusterFractionI + res.inClusterFractionV) / 2.0 << ",\n";
  outfile << "\"coords\": [";
  printDefects(res.defects);
  outfile << "],\n";
  outfile << "\"clusters\": ";
  outfile << "{";
  printClusterIds(res.clusters, outfile);
  outfile << "}";
  outfile << ",\n";
  outfile << "\"clusterSizes\": ";
  outfile << "{";
  printClusterIVs(res.clustersIV, outfile);
  outfile << "}";
  outfile << ",\n";
  outfile << "\"features\": ";
  outfile << "{";
  printFeats(res.feats, outfile);
  outfile << "}";
  outfile << ",\n";
  outfile << "\"coDefects\": [";
  avi::writeVector(res.coDefects, outfile);
  outfile << "],\n";
  outfile << "\"extraCols\": [";
  printExtraCols(res.extraCols, outfile);
  outfile << "]\n";
}

void avi::infoToKeyValue(std::ostream &outfile,
                              const avi::InputInfo &i,
                              const avi::ExtraInfo &ei) {
  outfile << "\"xyzFilePath\": \"" << i.xyzFilePath << "\",\n"
          << "\"id\": \"" << ei.id << "\",\n"
          << "\"substrate\": \"" << ei.substrate << "\",\n"
          << "\"simulationCode\": \"" << strSimulationCode(i.xyzFileType)
          << "\",\n"
          << "\"energy\":" << ei.energy << ",\n"
          << "\"structure\": \"" << i.structure << "\",\n"
          << "\"simulationTime\":" << ei.simulationTime << ",\n"
          << "\"ncell\":" << i.ncell << ",\n"
          << "\"boxSize\":" << i.boxSize << ",\n"
          << "\"origin\": [" << i.originX << ", " << i.originY << ", "
          << i.originZ << "],\n"
          << "\"rectheta\":" << ei.rectheta << ",\n"
          << "\"recphi\":" << ei.recphi << ",\n"
          << "\"pka\": [" << ei.xrec << ", " << ei.yrec << ", "
          << ei.zrec << "],\n"
          << "\"xyzColumnStart\":" << i.xyzColumnStart << ",\n"
          << "\"extraColumns\": [" << i.extraColumnStart << ", " << i.extraColumnEnd << "],\n"
          << "\"latticeConst\":" << i.latticeConst << ",\n"
          << "\"temperature\":" << i.temperature << ",\n"
          << "\"infile\": \"" << ei.infile << "\",\n"
          << "\"tags\": \"" << ei.tags << "\",\n"
          << "\"potentialUsed\": \"" << ei.potentialUsed << "\",\n"
          << "\"author\": \"" << ei.author << "\",\n"
          << "\"isPkaGiven\":" << ei.isPkaGiven << ",\n"
          << "\"es\":" << ei.es << ",\n"
          << "\"originType\":" << i.originType; // << ",\n";
}

void avi::configToKeyValue(std::ostream &outfile,
                                const avi::Config &c) {
  outfile << "\"version\": \"" << "0.4" << "\",\n"
          << "\"onlyDefects\": \"" << c.onlyDefects << "\",\n"
          << "\"isFindDistribution\": \"" << c.isFindDistribAroundPKA << "\",\n"
          << "\"isFindClusterFeatures\": \"" << c.isFindClusterFeatures << "\",\n"
          << "\"isIgnoreBoundaryDefects\": \"" << c.isIgnoreBoundaryDefects << "\",\n"
          << "\"isAddThresholdDefects\": \""
          << c.isAddThresholdDefects << "\",\n"
          << "\"filterZeroSizeClusters\":" << c.filterZeroSizeClusters << ",\n"
          << "\"thresholdFactor\":" << c.thresholdFactor << ",\n"
          << "\"logMode\":" << c.logMode << ",\n"
          << "\"logFilPath\": \"" << c.logFilePath << "\",\n"
          << "\"outputJSONFilePath\": \"" << c.outputJSONFilePath << "\""; // << ",\n";
}

void avi::printJson(std::ostream &outfile, const avi::InputInfo &i,
                         const avi::ExtraInfo &ei,
                         const avi::resultsT &res) {
  outfile << std::fixed << std::setprecision(4);
  outfile << "{";
  avi::infoToKeyValue(outfile, i, ei);
  outfile << ",\n";
  avi::resToKeyValue(outfile, res);
  outfile << "}\n";
}
// const char * c = outfile.str().c_str();
