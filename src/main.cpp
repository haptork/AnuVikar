#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <clipp.h>

#include <cluster2features.hpp>
#include <helper.hpp>
#include <logger.hpp>
#include <printJson.hpp>
#include <reader.hpp>
#include <infoReader.hpp>

/*
Reads command line arguments using clipp
*/
auto getConfig(int argc, char *argv[]) {
  using clipp::option;
  using clipp::value;
  avi::Config config;
  std::vector<std::string> files;
  int &logMode = config.logMode;
  auto help = false;
  std::string name;
  auto cli =
      (option("-sm")
           .set(config.onlyDefects, true)
           .doc("Switch: Compute minimum info: only the defect coordinates"),
       option("-sf")
           .set(config.isFindClusterFeatures, false)
           .doc("Switch: Do not add cluster features for patthern matching & "
                "classification."),
       option("-sb")
           .set(config.isIgnoreBoundaryDefects, false)
           .doc("Switch: Do not ignore boundary defects. Useful if defects can appear "
                "at boundary due to PBC if origin / offset is not given "
                "properly in MD simulations."),
       option("-sz")
           .set(config.filterZeroSizeClusters, true)
           .doc("Switch: Ignore clusters that appear due to purely threshold "
                "based interstitial-vacancy pairs and can annihilate totally."),
       option("-sd")
           .set(config.isAddThresholdDefects, false)
           .doc("Switch: Do not add threshold based displaced atom-lattice site pairs."),
       option("-sc")
           .set(config.safeRunChecks, false)
           .doc("Switch: Do not check and ignore files with anomalous number "
                "or proportion of defects."),
       option("-pf") &
           value("frames", config.allFrames)
           .doc("Switch: If multiple frames to be processed, total no. of frames. (default is 0 i.e. last frame only)"),
       option("-pt") &
           value("threshold", config.thresholdFactor)
               .doc("param (default 0.345): threshold factor for threshold "
                    "based interstitials (threshold value will be factor * "
                    "latticConst) if not switched off with -st"),
       option("-pc") &
           value("safety factor", config.thresholdFactor)
               .doc("param (default 50): safety factor for checks, lower value "
                    "implies stricter checks to ignore files. Only matters if "
                    "safety checks are not disabled by using -sc."),
       option("-o") & value("out file", config.outputJSONFilePath)
                          .doc("(default cascades-data.json), output JSON file "
                               "name / path."),
       option("-lf") &
           value("log file", config.logFilePath)
               .doc("log (default log-anuvikar-cpp.txt): file name"),
       option("-ld")
           .set(logMode, logMode | avi::LogMode::debug)
           .doc("log: Enable logging for debug"),
       option("-li")
           .set(logMode, logMode | avi::LogMode::info)
           .doc("log: Enable logging for info"),
       option("-lw")
           .set(logMode, logMode | avi::LogMode::warning)
           .doc("log: Enable logging for warning"),
       option("-le")
           .set(logMode, logMode | avi::LogMode::error)
           .doc("log: Enable logging for error"),
       option("-ln")
           .set(logMode, avi::LogMode::none | avi::LogMode::none)
           .doc("log: Disable logging"),
       option("-la")
           .set(logMode, avi::LogMode::all | avi::LogMode::all)
           .doc("log: enable all "),
       clipp::any_other(files), option("-help").set(help));
  if (clipp::parse(argc, argv, cli) && !help && !files.empty()) {
    return std::make_tuple(config, files, true);
  }
  std::cout << clipp::usage_lines(cli, "Csaransh") << '\n';
  if (files.empty()) {
    std::cout << "No input xyz file provided as cmd argument.\n";
  } else {
    std::cout << "xyz files as other cmd argument.\n";
  }
  std::cout << clipp::documentation(cli) << '\n';
  return std::make_tuple(config, files, false);
}

/*
Cmd-line args are parsed to get config. Each xyz file is processed and 
results are appended to the output file.
*/
int main(int argc, char *argv[]) {
  using avi::Logger; // Logger::inst().mode(avi::LogMode::warning | avi::LogMode::error);
  avi::Config config;
  std::vector<std::string> files;
  bool isConfigParse;
  // parsing cmd-line arguments to get config parameters
  std::tie(config, files, isConfigParse) =
      getConfig(argc, argv); // TODO: cook it using cmd-line flags
  if (!isConfigParse) return 1;
  // Using config for log and output-file
  Logger::inst().mode(config.logMode);
  Logger::inst().file(config.logFilePath);
  const std::string outpath{config.outputJSONFilePath};
  std::ofstream outfile{outpath};
  if (!outfile.is_open()) {
    std::cerr << "The output path " + outpath + " is not accessible.\n";
    return 1;
  }
  outfile << "{\"meta\": {\n";
  avi::configToKeyValue(outfile, config);
  outfile << "}\n, \"data\": [\n";
  std::cout << "Total files to process: " << (files.size()) << '\n'
            << std::flush;
  Logger::inst().log_info("Total files to process: " +
                          std::to_string(files.size()));
  Logger::inst().log_info("Started writing to output file \"" + outpath + "\"");
  auto success = 0;
  int curIndex = 0;
  avi::InputInfo info;
  avi::ExtraInfo extraInfo;
  auto infoStatus = avi::infoFrom::stdinOrFile;
  for (const auto &file : files) {
    std::cout << "\rCurrently processing file " << curIndex + 1 << std::flush;
    Logger::inst().log_info("Started processing file \"" + file + "\"");
    avi::ErrorStatus ret;
    int curSuccess = 0;
    std::tie(ret, curSuccess) = processFileTimeCmd(file, outfile, config, success, info, extraInfo, infoStatus);
    if (avi::ErrorStatus::noError != ret) {
      std::cerr << "\nError in processing file " << file << '\n';
      std::cerr << errToStr(ret) << '\n';
      Logger::inst().log_error("Error in processing file \"" + file + "\" " +
                               errToStr(ret));
    } else {
      if (config.allFrames) Logger::inst().log_info("Finished processing" + std::to_string(curSuccess) +" frames in file \"" + file + "\"");
      success += curSuccess;
      if (curIndex != files.size() - 1) outfile << ",";
      outfile << "\n";
    }
    curIndex++;
  }
  outfile << "]}"
          << "\n";
  outfile.close();
  std::cout << '\r' << curIndex << " out of " << curIndex
            << " processed successfully.\n";
  std::cout << "Output file written " + outpath << '\n';
  Logger::inst().log_info("Output file written " + outpath);
  return 0;
}
