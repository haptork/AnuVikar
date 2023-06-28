#include <algorithm>
#include <numeric>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <AddOffset.hpp>
#include <NextExpected.hpp>
#include <helper.hpp>
#include <logger.hpp>
#include <optimizeOrigin.hpp>
#include <results.hpp>
#include <xyz2defects.hpp>
#include <xyzReader.hpp>

#include <iostream>
#include <unordered_set>

// TODO import from printJson
auto strSimulationCodeD(avi::XyzFileType code) {
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



auto getThresh(const avi::InputInfo &info, const double &factor) {
  auto tempFactor = 1.0 + (info.temperature - 900)/10000;
  if (tempFactor < 1.0) tempFactor = 1.0;
  return (factor * info.latticeConst) * tempFactor ;
}

// tells if two coordinates are equal within epsilon = 1e-6 range
auto cmpApprox(const std::array<double, 3> &c1,
               const std::array<double, 3> &c2) {
  // constexpr auto epsilon = std::numeric_limits<double>::epsilon();
  auto epsilon = 1e-2;
  for (auto i : {0, 1, 2}) {
    if (fabs(c1[i] - c2[i]) > epsilon) return (c1[i] > c2[i]) ? 1 : -1;
  }
  return 0;
}

// tells if a coordinate exists in vector of coordinates. The equality is
// measured by cmpApprox
auto existAnchorVac(const std::array<double, 3> &t,
                    const std::vector<std::array<double, 3>> &v) {
  auto it = std::lower_bound(
      begin(v), end(v), t,
      [](const std::array<double, 3> &a, const std::array<double, 3> &b) {
        if (cmpApprox(a, b) < 0) return true;
        return false;
      });
  if (it != std::end(v) && cmpApprox(*it, t) == 0) return true;
  return false;
}

// tells if a coordinate exists in vector of coordinates + bool. The equality is
// measured by cmpApprox
auto existAnchorInter(
    const std::array<double, 3> &t,
    const std::vector<std::tuple<std::array<double, 3>, bool>> &v) {
  auto temp = std::make_tuple(t, false);
  auto it =
      std::lower_bound(begin(v), end(v), temp,
                       [](const std::tuple<std::array<double, 3>, bool> &a,
                          const std::tuple<std::array<double, 3>, bool> &b) {
                         if (cmpApprox(std::get<0>(a), std::get<0>(b)) < 0)
                           return true;
                         return false;
                       });
  if (it != std::end(v) && cmpApprox(std::get<0>(*it), t) == 0) return true;
  return false;
}

// sets the surviving flag of an interstitial and vacancy if they are closer
// than a threshold
auto clean(
    std::vector<std::tuple<avi::Coords, avi::Coords, bool, size_t>> &inter,
    std::vector<std::tuple<avi::Coords, bool>> &vac, double latticeConst) {
  const auto &thresh = latticeConst; // * sqrt(3) / 2;
  for (size_t i = 0; i < vac.size(); ++i) {
    if (!std::get<1>(vac[i])) continue;
    auto min = thresh + 1e-6;
    size_t minj = 0;
    for (size_t j = 0; j < inter.size(); ++j) {
      if (!std::get<2>(inter[j])) continue;
      auto dist =
          avi::calcDist(std::get<0>(vac[i]), std::get<1>(inter[j]));
      if (dist < min) {
        min = dist;
        minj = j;
      }
    }
    if (min < thresh) {
      std::get<1>(vac[i]) = false;
      std::get<2>(inter[minj]) = false;
    }
  }
}

std::tuple<avi::xyzFileStatus, avi::dispCoords, std::vector<std::vector<double>>>
getDisplacedAtomsTime(avi::InputInfo &info, avi::ExtraInfo &extraInfo,
         const avi::Config &config, std::istream &infile, avi::frameStatus &fs) {
  using avi::Coords;
  using std::get;
  using std::string;
  using std::tuple;
  using std::vector;
  std::tuple<avi::xyzFileStatus, avi::dispCoords, std::vector<std::vector<double>>> res;
  auto &atoms = std::get<1>(res);
  //std::array<std::vector<avi::Coords>, 2> atoms;
  // const auto latConst = info.latticeConst;
  std::string line;
  // read file and apply object
  std::array<avi::Coords, 2> c;
  std::vector<double> ec;
  avi::lineStatus ls;
  auto origin = std::array<double,3>{{info.originX, info.originY, info.originZ}};
  auto latConst = info.latticeConst;
  auto obj = avi::AddOffset{latConst, info.structure, origin};
  auto fn = [&latConst] (double x) { return x * latConst; };
  std::get<0>(res) = avi::xyzFileStatus::eof;
  while (std::getline(infile, line)) {
    std::tie(ls, c, ec) = avi::getCoordDisplaced(line, info.extraColumnStart, info.extraColumnEnd);
    if (ls == avi::lineStatus::coords &&
        fs == avi::frameStatus::inFrame) {
      auto vacC = std::get<0>(obj(c[0]));
      std::transform(begin(vacC), end(vacC), begin(vacC), fn);
      atoms.vacs.emplace_back(std::move(vacC));
      atoms.sias.emplace_back(std::make_pair(c[1], atoms.sias.size()));
      if(!ec.empty())  get<2>(res).emplace_back(ec);
    } else if (ls == avi::lineStatus::frameBorder) {
      if (fs != avi::frameStatus::prelude) {
        fs = avi::frameStatus::inFrame;
        if (config.allFrames && !atoms.vacs.empty()) {
          std::get<0>(res) = avi::xyzFileStatus::reading;
          break;
        }
        atoms.vacs.clear();
        atoms.sias.clear();
        get<2>(res).clear();
        avi::Logger::inst().log_warning(info.xyzFilePath + ", " + extraInfo.infile + ": " + " Multiple frames in file. Reading last one.");
      }
      fs = avi::frameStatus::inFrame;
    }
  }
  return res;
}

std::tuple<avi::xyzFileStatus, std::vector<avi::offsetCoords>, std::vector<std::vector<double>>>
getAtomsTime(avi::InputInfo &info, avi::ExtraInfo &extraInfo,
         const avi::Config &config, std::istream &infile, avi::frameStatus &fs) {
  using avi::Coords;
  using avi::strAr;
  using std::get;
  using std::string;
  using std::tuple;
  using std::vector;
  vector<Coords>
      atoms; // for all the atoms along with nearest closest sites and offset
  std::tuple<avi::xyzFileStatus, vector<avi::offsetCoords>, std::vector<std::vector<double>>>
      res; // for all the atoms along with nearest closest sites and offset
  const auto &latConst = info.latticeConst;
  if (info.ncell > 0) atoms.reserve(info.ncell * info.ncell * info.ncell * 2);
  std::string line;
  // read file and apply object
  Coords c;
  std::vector<double> ec;
  avi::lineStatus ls;
  std::get<0>(res) = avi::xyzFileStatus::eof;
  //std::cout<<"xyzColumn start " << info.xyzColumnStart << '\n';
  //std::cout<<"extraColumn start " << info.extraColumnStart << '\n';
  while (std::getline(infile, line)) {
    std::tie(ls, c, ec) = avi::getCoord(line, fs, info, extraInfo);
    // TODO: add ec to props
    if (ls == avi::lineStatus::coords &&
        fs == avi::frameStatus::inFrame) {
      atoms.emplace_back(c);
      get<2>(res).emplace_back(ec);
    } else if (ls == avi::lineStatus::inFrameCoords) {
      atoms.emplace_back(c);
      get<2>(res).emplace_back(ec);
    } else if (ls == avi::lineStatus::frameBorder) {
      if (fs != avi::frameStatus::prelude) {
        fs = avi::frameStatus::inFrame;
        if (config.allFrames && !atoms.empty()) {
          std::get<0>(res) = avi::xyzFileStatus::reading;
          break;
        }
        if (!atoms.empty()) {
          //std::cout << "Atoms before clearing: " << atoms[0][0] << ", " << atoms[0][1] << ", " << atoms[0][2] << '\n';
          avi::Logger::inst().log_warning(info.xyzFilePath + ", " + extraInfo.infile + ": " + " Multiple frames in file. Reading last one.");
        }
        atoms.clear();
      }
      fs = avi::frameStatus::inFrame;
        //atoms.clear(); // Ignoring the frame before this one
    }
  }
  //std::cout <<"atoms here: " << atoms.size() << '\n';
  if (atoms.empty()) return res;
  //std::cout <<"first atom: " << atoms[0][0] << ", " << atoms[0][1] << ", " << atoms[0][2] << '\n';
  // assuming bcc /fcc structure and perfect initial. TODO: for imperfect skip
  if (info.structure[0] == 'b') info.ncell = std::round(std::cbrt(atoms.size() / 2.0));
  else if (info.structure[0] == 'f') info.ncell = std::round(std::cbrt(atoms.size() / 4.0));
  auto secLatConst = (info.boxSize > 0.0) ? info.boxSize / info.ncell : -1.0;
  if (secLatConst > 0.0 && std::fabs(info.latticeConst - secLatConst) < 1e-3) {
    secLatConst = -1.0;
  }
  if (info.latticeConst < 0.0) {
    info.latticeConst = secLatConst;
    secLatConst = -1.0;
  }
  //std::cout << "lat const: " << info.latticeConst << '\n';
  //std::cout << "origin type : " << info.originType << '\n';
  std::vector<std::tuple<avi::Coords, double, std::string, int>> combos;
  if (info.originType > 0) {
    auto originEstimated = avi::estimateOrigin(atoms, info.latticeConst);
    combos.emplace_back(
        avi::Coords{
            {originEstimated[0], originEstimated[1], originEstimated[2]}},
        info.latticeConst, "estimated origin, given latConst", 0);
    if (secLatConst > 0.0) {
      auto originEstimatedSec = avi::estimateOrigin(atoms, secLatConst);
      combos.emplace_back(
          avi::Coords{{originEstimatedSec[0], originEstimatedSec[1],
                            originEstimatedSec[2]}},
          secLatConst, "estimated originSec, boxSize/ncell latconst", 0);
    }
  }
  if ((secLatConst > 0.0 && info.originType == 0) || info.originType == 2) {
    combos.emplace_back(
        avi::Coords{{info.originX, info.originY, info.originZ}},
        info.latticeConst, "given origin, given latconst", 0);
    if (secLatConst > 0.0)
      combos.emplace_back(
          avi::Coords{{info.originX, info.originY, info.originZ}},
          secLatConst, "given origin, boxSize/ncell latconst", 0);
  }
  if (!combos.empty()) {
    auto leastIndex = 0;
    if (combos.size() > 1) {
      auto leastInterThresh = std::numeric_limits<int>::max();
      auto thresh = getThresh(info, config.thresholdFactor);
      for (size_t i = 0; i < combos.size(); i++) {
        auto curOrigin = std::get<0>(combos[i]);
        auto curLatConst = std::get<1>(combos[i]);
        auto obj = avi::AddOffset{curLatConst, info.structure, curOrigin};
        std::get<3>(combos[i]) = std::count_if(
            begin(atoms), end(atoms), [&obj, thresh](const auto &it) {
              return std::get<1>(obj(it)) > thresh;
            });
        avi::Logger::inst().log_debug(
            std::get<2>(combos[i]) + " origin, latconst " + strAr(curOrigin) +
            ", " + std::to_string(curLatConst) + " : nThreshold " +
            std::to_string(std::get<3>(combos[i])));
        if (std::get<3>(combos[i]) < leastInterThresh) {
          leastIndex = i;
          leastInterThresh = std::get<3>(combos[i]);
        }
      }
    }
    if (combos.size() > 1 && info.originType == 1 &&  // estimated origin but two lattice constants, one given, one from boxdimensions & unit-cells
         (secLatConst > 0.0 && std::fabs(info.latticeConst - secLatConst) > 1e-2)) {
      auto mx = std::max(std::get<3>(combos[0]), std::get<3>(combos[1]));
      if (mx > 1.1 * std::get<3>(combos[leastIndex])) {
        const std::string msgPre = (leastIndex == 0) ? "given value of box-size might not be exact." : "given value of lattice constant might not be exact";
        if (leastIndex == 0) {
          avi::Logger::inst().log_warning(info.xyzFilePath + ", " + extraInfo.infile + ": " + msgPre + "Difference in given lattice constant & box-size based lattice constant caused difference in"
               "displaced atom values: " + std::to_string(std::get<3>(combos[0])) + ", " + std::to_string(std::get<3>(combos[1])) +
               "lattice constants (given, calc): " + std::to_string(info.latticeConst) + ", " + std::to_string(secLatConst));
        } else {
          avi::Logger::inst().log_warning(info.xyzFilePath + ", " + extraInfo.infile + ": reverse " + msgPre + "Difference in given lattice constant & box-size based lattice constant caused difference in"
               "displaced atom values: " + std::to_string(std::get<3>(combos[0])) + ", " + std::to_string(std::get<3>(combos[1])) +
               "lattice constants (given, calc): " + std::to_string(info.latticeConst) + ", " + std::to_string(secLatConst));
        }
      }
    }
    info.originX = std::get<0>(combos[leastIndex])[0];
    info.originY = std::get<0>(combos[leastIndex])[1];
    info.originZ = std::get<0>(combos[leastIndex])[2];
    info.latticeConst = std::get<1>(combos[leastIndex]);
  }
  std::get<1>(res).reserve(atoms.size());
  auto origin = avi::Coords{{info.originX, info.originY, info.originZ}};
  auto obj = avi::AddOffset{info.latticeConst, info.structure, origin};
  //std::cout << "latInfo: " << info.latticeConst << ", " << info.structure << ", " << origin[0] << ", " << origin[1] << ", " << origin[2] << '\n';
  std::transform(begin(atoms), end(atoms), std::back_inserter(std::get<1>(res)), obj);
  //std::cout << "\natoms1: " << atoms[0][0] << " | " << atoms[atoms.size() - 1][0] << '\n';
  //std::cout << "\natoms res: " << std::get<3>(std::get<1>(res)[0]) << " | " << std::get<3>(std::get<1>(res)[1]) << "\n";
  if (info.boxSize < 0.0) { info.boxSize = info.latticeConst * info.ncell; }
  return res;
}

avi::DefectRes
avi::xyz2defectsTime(avi::InputInfo &mainInfo,
                      avi::ExtraInfo &extraInfo,
                      const avi::Config &config, std::istream &infile, avi::frameStatus &fs) {
  auto atoms = getAtomsTime(mainInfo, extraInfo, config, infile, fs);
  if (std::get<1>(atoms).empty() && mainInfo.xyzFileType != avi::XyzFileType::generic) {
    infile.clear();
    infile.seekg(0);
    auto temp = mainInfo.xyzFileType;
    mainInfo.xyzFileType = avi::XyzFileType::generic;
    if (mainInfo.xyzColumnStart < 0) mainInfo.xyzColumnStart = 0;
    atoms = getAtomsTime(mainInfo, extraInfo, config, infile, fs);
    if (!std::get<1>(atoms).empty()) {
      Logger::inst().log_warning(extraInfo.infile +": " + mainInfo.xyzFilePath + ": Default file format " + strSimulationCodeD(temp) + " is not correct. Fallback to generic reading works.");
    }
  }
  //std::cout << "\natoms: " << atoms.second.size() << '\n';
  //std::cout << "\nxyzCol: " << mainInfo.xyzColumnStart << '\n';
  if (std::get<1>(atoms).empty())
    return std::make_tuple(std::get<0>(atoms), ErrorStatus::xyzFileReadError, avi::DefectVecT{}, std::vector<int>{}, std::vector<std::vector<double>>{});
  return  atoms2defects(atoms, mainInfo, extraInfo, config, (mainInfo.structure == "bcc"));
}

avi::DefectRes
avi::displaced2defectsTime(avi::InputInfo &mainInfo,
                      avi::ExtraInfo &extraInfo,
                      const avi::Config &config, std::istream &infile, avi::frameStatus &fs) {
  auto atoms = getDisplacedAtomsTime(mainInfo, extraInfo, config, infile, fs);
  if (std::get<1>(atoms).vacs.empty())
    return std::make_tuple(std::get<0>(atoms), ErrorStatus::noError, avi::DefectVecT{}, std::vector<int>{}, std::vector<std::vector<double>>{});
  return avi::displacedAtoms2defects(atoms, mainInfo.latticeConst, mainInfo.boxSize);
}

/*
void testIt(std::vector<std::tuple<avi::Coords, double, avi::Coords>>
atoms, avi::NextExpected n) { std::ofstream f{"junk_got.txt"}; for (auto it
: atoms) { for (auto jt: std::get<0>(it)) { f << jt << " ";
    }
    f << " | ";
    for (auto jt: std::get<2>(it)) {
      f << jt << " ";
    }
    f << '\n';
  }
  f.close();
  std::ofstream fe{"junk_expected.txt"};
  for (auto jt: n.cur()) {
    fe << jt << " ";
  }
  fe << '\n';
  while (!n.allMax()) {
   n.increment();
   for (auto jt: n.cur()) {
      fe << jt << " ";
    }
    fe << " | ";
   for (auto jt: n.minCur()) {
      fe << jt << " ";
    }
    fe << " | ";
    for (auto jt: n.maxCur()) {
      fe << jt << " ";
    }
    fe << '\n';
   }
  fe.close();
}
*/

avi::Coords avi::getInitialMax(const avi::Coords &origin,
                                         const avi::Coords &maxes) {
  auto maxes1 = maxes;
  constexpr auto epsilon = 1e-2;
  for (auto i = 0; i < 3; i++) {
    if (maxes[i] - origin[i] - (int(maxes[i] - origin[i])) < epsilon) {
      maxes1[i] = maxes[i];
      // maxes2[i] = maxes[i] - 0.5;
    } else {
      maxes1[i] = maxes[i] - 0.5;
      // maxes2[i] = maxes[i];
    }
  }
  return maxes1;
}

avi::Coords avi::getInitialMaxFcc(const avi::Coords &origin,
                                         const avi::Coords &maxes) {
  auto maxes1 = maxes;
  constexpr auto epsilon = 1e-2;
  //for (auto i = 2; i < 3; i++) {
  auto i = 2;
    if (maxes[i] - origin[i] - (int(maxes[i] - origin[i])) < epsilon) {
      maxes1[i] = maxes[i];
      // maxes2[i] = maxes[i] - 0.5;
    } else {
      maxes1[i] = maxes[i] - 0.5;
      // maxes2[i] = maxes[i];
    }
  //}
  return maxes1;
}

void setCenter(avi::ExtraInfo &extraInfo, avi::Coords minC,
               avi::Coords maxC, double latConst) {
  if (!extraInfo.isPkaGiven) {
    extraInfo.xrec = ((maxC[0] - minC[0]) / 2.0 + minC[0])*latConst;
    extraInfo.yrec = ((maxC[1] - minC[1]) / 2.0 + minC[1])*latConst;
    extraInfo.zrec = ((maxC[2] - minC[2]) / 2.0 + minC[2])*latConst;
  }
}

avi::DefectRes avi::atoms2defects(
    std::tuple<avi::xyzFileStatus, std::vector<avi::offsetCoords>, std::vector<std::vector<double>>> stAtoms,
    avi::InputInfo &info, avi::ExtraInfo &extraInfo,
    const avi::Config &config, bool isBcc) {
  using avi::Coords;
  using std::get;
  using std::tuple;
  using std::vector;
  bool isExtraCol = (!(std::get<2>(stAtoms)).empty() && !std::get<2>(stAtoms)[0].empty());
  bool isIncludeId = isExtraCol || info.extraColumnStart == -2;
  auto &atoms = std::get<1>(stAtoms);
  std::sort(begin(atoms), end(atoms));
  // std::cout << "\natoms: " << std::get<0>(atoms[0])[0] << " | " << std::get<0>(atoms[atoms.size() - 1])[0] << '\n';
  const auto minmax = std::minmax_element(
      begin(atoms), end(atoms), [](const auto &ao, const auto &bo) {
        const auto &a = std::get<0>(ao);
        const auto &b = std::get<0>(bo);
        return (a[0] + a[1] + a[2]) < (b[0] + b[1] + b[2]);
      });
  auto firstRow = std::get<0>(*(minmax.first));
  auto lastRow = std::get<0>(*(minmax.second));
  auto maxes = lastRow;
  auto maxesInitial = (isBcc) ? avi::getInitialMax(firstRow, maxes)
                              : avi::getInitialMaxFcc(firstRow, maxes);
  auto maxesFinal = (isBcc) ? maxes : maxesInitial;
  avi::NextExpected nextExpected{firstRow, maxes, maxesInitial, maxesFinal};
  setCenter(extraInfo, firstRow, lastRow, info.latticeConst);
  std::tuple<double, Coords, Coords, size_t> pre;
  bool isPre = false;
  const auto latConst = info.latticeConst;
  const auto thresh = getThresh(
      info, config.thresholdFactor); // for rough interstitial threshold
  vector<tuple<double, Coords, Coords, size_t>>
      interThresh; // interstitials based on rough thresh
  vector<tuple<Coords, Coords, bool, size_t>> interstitials;
  vector<tuple<Coords, bool>> vacancies;
  std::vector<int> vacSias;
  std::vector<std::vector<double>> ids;
  vector<tuple<Coords, Coords, bool, size_t>> freeInterstitials;
  const auto boundaryThresh = (config.isIgnoreBoundaryDefects) ? 0.501 : 0.00;
  auto boundaryPred = [](Coords ar, Coords min, Coords max, double thresh) {
    for (size_t i = 0; i < ar.size(); i++) {
      if (ar[i] - thresh < min[i] || ar[i] + thresh > max[i]) { return true; }
    }
    return false;
  };
  const auto &extraFactor = config.extraDefectsSafetyFactor;
  std::string errMsg = "";
  if (!(Logger::inst().mode() & LogMode::info))
    errMsg = "from \"" + info.xyzFilePath + "\", \"" + extraInfo.infile + "\": ";
  for (const auto &row : atoms) {
    Coords ar = get<0>(row);
    double curOffset = get<1>(row);
    Coords c = std::get<2>(row);
    size_t id = std::get<3>(row);
    const auto cmpRes = cmpApprox(ar, nextExpected.cur());
    if (nextExpected.allMax()) break;
    // threshold
    // std::cout << ar[0] << ", " << config.isAddThresholdDefects << '\n';
    if (config.isAddThresholdDefects && curOffset > thresh) {
      interThresh.emplace_back(curOffset, ar, c, get<3>(row));
    }
    // clean
    using vecB = std::vector<std::tuple<Coords, double, Coords>>;
    if (cmpRes == 0) { // atom at correct lattice site
      if (isBcc) nextExpected.increment(); else nextExpected.incrementFcc();
      pre = std::make_tuple(curOffset, ar, c, id);
      isPre = true;
    } else if (cmpRes > 0) { // at site after expected
      vecB res;
      while (cmpApprox(ar, nextExpected.cur()) > 0 &&
             cmpApprox(nextExpected.cur(), lastRow) <= 0 &&
             !nextExpected.allMax()) {
        if (!boundaryPred(nextExpected.cur(), nextExpected.minCur(),
                          nextExpected.maxCur(), boundaryThresh)) {
          vacSias.emplace_back(interstitials.size());
          vacancies.emplace_back(nextExpected.cur(), true);
          //avi::Logger::inst().log_debug("vac: " + strAr(nextExpected.cur()) + " coord: " + strAr(ar)); 
          if (config.safeRunChecks &&
              vacancies.size() * extraFactor > atoms.size()) {
            Logger::inst().log_error( errMsg + 
                "not processing since too many vacancies, this might be due to "
                "wrong inputs like latticeConst, boxDim or corrupt xyz file. "
                "(v): " +
                std::to_string(vacancies.size()));
            return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::vacOverflow, DefectVecT{}, vector<int>{}, ids);
          }
        }
        if (isBcc) nextExpected.increment(); else nextExpected.incrementFcc();
      }
      if (isBcc) nextExpected.increment(); else nextExpected.incrementFcc();
    } else { // atom sits at the same lattice site again
      vecB res;
      if (isPre && std::get<1>(pre) == ar) {
        auto isPreReal = std::get<0>(pre) > curOffset;
        vacancies.emplace_back(ar, false); // dummy vacancy added
        interstitials.emplace_back(std::get<1>(pre), std::get<2>(pre),
                                   isPreReal, std::get<3>(pre));
        interstitials.emplace_back(
            ar, c,
            !isPreReal, id); // both interstitials for structures like dumbbells
        vacSias.push_back(interstitials.size());
      } else {
        if (boundaryPred(ar, nextExpected.minCur(), nextExpected.maxCur(),
                         boundaryThresh) == false) {
          if (vacSias.size() > 0 && cmpApprox(std::get<1>(pre), ar) == 0) {
            interstitials.emplace_back(ar, c, true, id);
            vacSias[vacSias.size() - 1] = interstitials.size();
          } else {
            freeInterstitials.emplace_back(ar, c, true, id);
            // std::cout << "\n interstitial wo pre: " << c[0] << c[1] << c[2] << std::endl; 
          }
        }
      }
      if (config.safeRunChecks &&
          interstitials.size() * extraFactor > atoms.size()) {
        Logger::inst().log_error( errMsg + 
            "not processing since too many interstitials, this might be due to "
            "wrong inputs like latticeConst, boxDim or corrupt xyz file. "
            "(i): " +
            std::to_string(interstitials.size()));
        return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::siaOverflow, DefectVecT{}, vector<int>{}, ids);
      }
      isPre = false;
    }
  }
  if (config.safeRunChecks) {
    if (interstitials.size() * extraFactor > atoms.size() ||
        vacancies.size() * extraFactor > atoms.size()) {
      Logger::inst().log_error(errMsg +
          "not processing since too many interstitials and vacancies, this may "
          "be due to wrong inputs like latticeConst, boxDim or corrupt xyz "
          "file. (i, v): " +
          std::to_string(interstitials.size()) + ", " +
          std::to_string(vacancies.size()));
      return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::defectOverflow, DefectVecT{}, vector<int>{}, ids);
    }
    if (interThresh.size() > extraFactor * interstitials.size()) {
      Logger::inst().log_error(errMsg + 
          "not processing since number of threshold based displaced atoms are excessive, this may be due to "
          "inputs like latticeConst, boxDim or corrupt xyz file.  (i, v, "
          "interThresh): " +
          std::to_string(interstitials.size()) + ", " +
          std::to_string(vacancies.size()) + ", " +
          std::to_string(interThresh.size()));
      return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::threshOverflow, DefectVecT{}, vector<int>{}, ids);
    }
    if (interstitials.size() != vacancies.size()) {
      if ((int(interstitials.size()) / int(vacancies.size())) > extraFactor) {

        avi::Logger::inst().log_warning(errMsg +
            "sia and vacancy counts are different (i, v, "
            "interTresh): " +
            std::to_string(interstitials.size()) + ", " +
            std::to_string(vacancies.size()) + ", " +
            std::to_string(interThresh.size()));
      } else {
        avi::Logger::inst().log_debug(
            "sia and vacancy counts are different (i, v, "
            "interThresh): " +
            std::to_string(interstitials.size()) + ", " +
            std::to_string(vacancies.size()) + ", " +
            std::to_string(interThresh.size()));
      }
      if (interstitials.size() * (extraFactor / 10) < vacancies.size() ||
          vacancies.size() * (extraFactor / 10) < interstitials.size()) {

        avi::Logger::inst().log_error(errMsg+
            "not processing since the difference in sia and vacancy is too big (i, v, "
            "interThresh): " +
            std::to_string(interstitials.size()) + ", " +
            std::to_string(vacancies.size()) + ", " +
            std::to_string(interThresh.size()));
        return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::siaVacDiffOverflow, DefectVecT{}, vector<int>{}, ids);
      }
    }
  }
  // std::cout<<"\ninterstitials: " << interstitials.size() << '\n';
  // std::cout<<"\nvacancies: " << vacancies.size() << '\n';
  // std::cout<<"\ninterThresh: " << interThresh.size() << '\n';
  atoms.clear();
  avi::DefectVecT defects;
  defects.reserve(2 * vacancies.size());
  auto vacCleanSave = vacancies;
  auto anchor2coord = [&latConst](Coords c) {
    for (auto &x : c)
      x *= latConst;
    return c;
  };
  auto count = 1;
  for (auto &it : vacancies) {
    for (auto &x : std::get<0>(it))
      x *= latConst;
  }
  /*
  std::cout <<"\nis: \n";
  for(auto &it: interstitials) {
    std::cout<<get<1>(it)[0] << ", " << get<1>(it)[1] << '\n';
  }
  std::cout <<"\nvs: \n";
  for(auto &it: vacancies) {
    std::cout<<get<0>(it)[0] << ", " << get<0>(it)[1] << '\n';
  }
  std::cout <<"\nvacSias: \n";
  for(auto &it: vacSias) {
    std::cout<<it <<'\n';
  }
  */
  for (auto it: freeInterstitials) {
    interstitials.push_back(it);
  }
  clean(interstitials, vacancies, info.latticeConst); // clean again
  std::vector<int> vacSiasNu;
  for (auto i = 0; i < vacancies.size(); i++) {
    auto &it = vacancies[i];
    defects.emplace_back(std::move(std::get<0>(it)), false, count++, std::get<1>(it));
    //std::cout << "I is " << i << ": ";
    for (auto j = (i > 0) ? vacSias[i-1] : 0; j < vacSias[i]; j++) {
      //std::cout << j << ", ";
      auto &jt = interstitials[j];
      defects.emplace_back(std::move(get<1>(jt)), true, count++, get<2>(jt));
      if (isIncludeId) ids.push_back(std::vector<double>{double(get<3>(jt))});
      if (isExtraCol) {
        auto &dst = ids[ids.size()-1];
        auto &src = std::get<2>(stAtoms)[get<3>(jt)];
        dst.insert(dst.end(), src.begin(), src.end());
      }
    }
    //std::cout << '\n';
    vacSiasNu.push_back(defects.size());
  }
  /*
  std::cout <<"\ndefects: \n";
  for(auto &it: defects) {
    std::cout<<get<0>(it)[0] << ", " << get<0>(it)[1] << '\n';
  }
  */
  int extraDefects = 0;
  int maxExtraDefects = interstitials.size() * 2;
  if (maxExtraDefects > interThresh.size())
    maxExtraDefects = interThresh.size();
  std::nth_element(begin(interThresh), begin(interThresh) + maxExtraDefects,
                   end(interThresh), [](const auto &a, const auto &b) {
                     return std::get<0>(a) > std::get<0>(b);
                   });
  for (const auto &it : interThresh) {
    if (!existAnchorInter(std::get<1>(it), vacCleanSave)) {
      extraDefects++;
      if (config.safeRunChecks && extraDefects > maxExtraDefects) continue;
      defects.emplace_back(anchor2coord(std::move(get<1>(it))), false, count++,
                           false);
      defects.emplace_back(std::move(get<2>(it)), true, count++, false);
      if (isIncludeId) ids.push_back(std::vector<double>{double(get<3>(it))});
      if (isExtraCol) {
        auto &dst = ids[ids.size()-1];
        auto &src = std::get<2>(stAtoms)[get<3>(it)];
        dst.insert(dst.end(), src.begin(), src.end());
      }
      vacSiasNu.emplace_back(defects.size());
    }
  }
  for (const auto &jt : freeInterstitials) {
    defects.emplace_back(std::move(get<1>(jt)), true, count++, get<2>(jt));
    if (isIncludeId) ids.push_back(std::vector<double>{double(get<3>(jt))});
    if (isExtraCol) {
      auto &dst = ids[ids.size()-1];
      auto &src = std::get<2>(stAtoms)[get<3>(jt)];
      dst.insert(dst.end(), src.begin(), src.end());
    }
  }
  if (config.safeRunChecks && extraDefects > interstitials.size()) {

    if (extraDefects > maxExtraDefects) {
      Logger::inst().log_warning(errMsg+
          "Some threshold based defects ignored. Too many threshold based "
          "defects. (i, v, extraI): " +
          std::to_string(interstitials.size()) + ", " +
          std::to_string(vacancies.size()) + ", " +
          std::to_string(extraDefects));
    } else {
      Logger::inst().log_debug(
          "Too many threshold based defects. (i, v, extraI): " +
          std::to_string(interstitials.size()) + ", " +
          std::to_string(vacancies.size()) + ", " +
          std::to_string(extraDefects));
    }
  }
  return std::make_tuple(std::get<0>(stAtoms), ErrorStatus::noError, std::move(defects), std::move(vacSiasNu), std::move(ids));
}

auto cleanDisplaced(avi::DefectVecT &inter, avi::DefectVecT &vac,
                    double latticeConst) {
  using avi::DefectTWrap::coords;
  using avi::DefectTWrap::isSurviving;
  auto thresh = latticeConst; // * sqrt(3) / 2;
  for (size_t i = 0; i < vac.size(); ++i) {
    auto min = thresh + 1e-6;
    size_t minj = 0;
    for (size_t j = 0; j < inter.size(); ++j) {
      if (!isSurviving(inter[j])) continue;
      auto dist = avi::calcDist(coords(vac[i]), coords(inter[j]));
      if (dist < min) {
        min = dist;
        minj = j;
      }
    }
    if (min < thresh) {
      isSurviving(vac[i], false);
      isSurviving(inter[minj], false);
    }
  }
}

avi::DefectRes avi::displacedAtoms2defects(
    std::tuple<avi::xyzFileStatus, dispCoords, std::vector<std::vector<double>>> statoms, double latticeConst, double box) {
  using avi::Coords;
  using std::get;
  using std::tuple;
  using std::vector;
  bool isExtraCol = (!(std::get<2>(statoms)).empty() && !std::get<2>(statoms)[0].empty());
  bool isIncludeId = isExtraCol;// || info.extraColumnStart == -2;
  avi::DefectVecT defects;
  auto &atoms = std::get<1>(statoms);
  constexpr auto epsilon = 1e-4;
  const auto nn = (std::sqrt(3) * latticeConst) / 2 + epsilon;
  const auto threshDist = 0.345 * latticeConst + 1e-4;
  const auto threshDistSqr = threshDist * threshDist;
  const auto recombDist = latticeConst*1.5;//latticeConst; //0.345 * latticeConst + 1e-4;
  const auto recombDistSqr = recombDist * recombDist;
  const auto vacGroupDist = nn;
  const auto vacGroupDistSqr = vacGroupDist * vacGroupDist;
  const auto maxSqrDist = (recombDistSqr > vacGroupDistSqr) ? recombDistSqr : vacGroupDistSqr;
  std::vector<int> freeIs;
  const auto sortHelper = [](const std::array<double, 3> &c1, const std::array<double, 3> &c2) {
    return cmpApprox(c1, c2) < 0;
  };
  sort( atoms.vacs.begin(), atoms.vacs.end(), sortHelper);
  const auto areEq = [](const std::array<double, 3> &c1, const std::array<double, 3> &c2) {
    return cmpApprox(c1, c2) == 0;
  };
  atoms.vacs.erase(unique(atoms.vacs.begin(), atoms.vacs.end(), areEq), atoms.vacs.end());
  
  const auto& siaIn = atoms.sias;
  const auto& vacIn = atoms.vacs;
  const auto& threshBig = recombDistSqr;

  std::vector<double> distVac2;
  std::vector<int> vacOcc_;
  for (const auto& sia : siaIn) {
      double minDistance = std::numeric_limits<double>::max();
      int minIndex = -1;
      for (int i = 0; i < vacIn.size(); ++i) {
          double distance = calcDistSqr(sia.first, vacIn[i], box);
          if (distance < minDistance) {
              minDistance = distance;
              minIndex = i;
          }
      }
      distVac2.push_back(minDistance);
      vacOcc_.push_back(minIndex);
  }
  // Calculate vacOcc, freeSias, vacOccUnique
  std::vector<int> vacOcc;
  std::vector<int> freeSias;
  std::vector<int> vacOccUnique;
  std::unordered_set<int> vacOccUniqueSet;
  for (int i = 0; i < distVac2.size(); ++i) {
      if (distVac2[i] < threshBig) {
          vacOcc.push_back(vacOcc[i]);
          vacOccUniqueSet.insert(vacOcc[i]);
      } else {
          freeSias.push_back(i);
      }
  }
  std::vector<int> vac(vacIn.size());
  std::vector<int> vacRecomb;
  std::vector<int> indexOfVacRecomb;
  std::vector<double> vacRecombInverted;
  std::vector<int> sortOrder(vacRecombInverted.size());
  std::vector<int> sia(vacRecombInverted.size());
  //std::vector<int> vacRecombNonUnique(vacRecombInverted.size());

  auto vacEnd = std::copy_if(vacIn.begin(), vacIn.end(), vac.begin(),
                              [&](const avi::Coords& vacVec) {
                                  // TODO: fix this
                                  //int index = &vacVec - vacIn.data();
                                  return 1;//vacOccUniqueSet.find(index) == vacOccUniqueSet.end();
                              });
  vac.resize(std::distance(vac.begin(), vacEnd));

  std::copy_if(vacOccUnique.begin(), vacOccUnique.end(), std::back_inserter(vacRecomb),
                [&](int index) { return std::count(vacOcc.begin(), vacOcc.end(), index) > 1; });

  std::copy_if(vacOcc.begin(), vacOcc.end(), std::back_inserter(indexOfVacRecomb),
                [&](int index) { return std::find(vacRecomb.begin(), vacRecomb.end(), index) != vacRecomb.end(); });
                
  std::transform(indexOfVacRecomb.begin(), indexOfVacRecomb.end(), std::back_inserter(vacRecombInverted),
                  [&](int index) { return vacOcc[index]; });  

  std::iota(sortOrder.begin(), sortOrder.end(), 0);
  std::sort(sortOrder.begin(), sortOrder.end(),
              [&](int i1, int i2) { return vacRecombInverted[i1] < vacRecombInverted[i2]; });

  std::transform(sortOrder.begin(), sortOrder.end(), sia.begin(),
                   [&](int i) { return indexOfVacRecomb[i]; });

  //std::transform(sortOrder.begin(), sortOrder.end(), vacRecombNonUnique.begin(),
                   //[&](int i) { return vacRecombInverted[i]; });

  // TODO: initialize with sizes that are known
  std::vector<int> latticeSiteGroups;
  std::vector<std::vector<double>> ids;
  //latticeSiteGroups.reserve(100);
  for (int i = 0; i < vacRecomb.size(); i++) {
      const auto &vc = vacIn[vacRecomb[i]];
      defects.emplace_back(vc, false, defects.size()+1, false);
      // TODO: generalize for di-interstitials and more?
      auto dist1 = calcDistSqr(vc, siaIn[sia[i*2]].first, box);
      auto dist2 = calcDistSqr(vc, siaIn[sia[i*2+1]].first, box);
      defects.emplace_back(siaIn[sia[i*2]], false, defects.size()+1, dist1>=dist2);
      defects.emplace_back(siaIn[sia[i*2+1]], false, defects.size()+1, dist1<dist2);
      if (isIncludeId) {
        ids.push_back(std::vector<double>{double(sia[i*2])});
        ids.push_back(std::vector<double>{double(sia[i*2+1])});
      }
      if (isExtraCol) {
        auto &dst = ids[ids.size()-1];
        auto &src = std::get<2>(statoms)[ids[ids.size()-1][0]];
        dst.insert(dst.end(), src.begin(), src.end());
      }
      latticeSiteGroups.push_back(defects.size());
  }
  for (auto it : freeIs) {
    defects.emplace_back(atoms.sias[it].first, true, defects.size() + 1, true); // interstitial
    if (isIncludeId) ids.push_back(std::vector<double>{double(it)});
    if (isExtraCol) {
      auto &dst = ids[ids.size()-1];
      auto &src = std::get<2>(statoms)[ids[ids.size()-1][0]];
      dst.insert(dst.end(), src.begin(), src.end());
    }
  }
  return std::make_tuple(std::get<0>(statoms), avi::ErrorStatus::noError, std::move(defects), std::move(latticeSiteGroups), ids);
}

avi::DefectRes displacedAtoms2defectsOld(
    std::tuple<avi::xyzFileStatus, avi::dispCoords, std::vector<std::vector<double>>> statoms, double latticeConst, double box) {
  using avi::Coords;
  using std::get;
  using std::tuple;
  using std::vector;
  bool isExtraCol = (!(std::get<2>(statoms)).empty() && !std::get<2>(statoms)[0].empty());
  bool isIncludeId = isExtraCol;// || info.extraColumnStart == -2;
  avi::DefectVecT inter, vac, defects;
  auto &atoms = std::get<1>(statoms);
  constexpr auto epsilon = 1e-4;
  const auto nn = (std::sqrt(3) * latticeConst) / 2 + epsilon;
  const auto threshDist = 0.345 * latticeConst + 1e-4;
  const auto threshDistSqr = threshDist * threshDist;
  const auto recombDist = latticeConst;//latticeConst; //0.345 * latticeConst + 1e-4;
  const auto recombDistSqr = recombDist * recombDist;
  const auto vacGroupDist = nn;
  const auto vacGroupDistSqr = vacGroupDist * vacGroupDist;
  const auto maxSqrDist = (recombDistSqr > vacGroupDistSqr) ? recombDistSqr : vacGroupDistSqr;
  std::vector<int> freeIs;
  const auto sortHelper = [](const std::array<double, 3> &c1, const std::array<double, 3> &c2) {
    return cmpApprox(c1, c2) < 0;
  };
  sort( atoms.vacs.begin(), atoms.vacs.end(), sortHelper);
  const auto areEq = [](const std::array<double, 3> &c1, const std::array<double, 3> &c2) {
    return cmpApprox(c1, c2) == 0;
  };
  atoms.vacs.erase(unique(atoms.vacs.begin(), atoms.vacs.end(), areEq), atoms.vacs.end());
  std::vector<std::vector<std::pair<double, int>>> vacSias(atoms.vacs.size());
  std::vector<std::vector<int>> siaVacsFull(atoms.sias.size());
  for (auto is = 0; is < atoms.sias.size(); is++) {
    auto minDist = maxSqrDist;
    auto minDistIndex = -1;
    for (auto iv = 0; iv < atoms.vacs.size(); iv++) {
      auto disp = avi::calcDistSqr(atoms.sias[is].first, atoms.vacs[iv], box);
      if (disp < maxSqrDist) siaVacsFull[is].push_back(iv);// = disp;
      if (disp < minDist) {
        minDist = disp;
        minDistIndex = iv;
      }
    }
    if (minDistIndex >=0) {
      vacSias[minDistIndex].push_back(std::make_pair(minDist, is));
    } else {
      freeIs.push_back(is);
    }
  }
  for (auto i = 0; i < vacSias.size(); i++) {
    std::sort(begin(vacSias[i]), end(vacSias[i]));
  }
  /*
  std::cout << "siaFull: \n";
  for (auto is = 0; is < siaVacsFull.size(); is++) {
    for (const auto &jt : atoms[1][is]) std::cout << jt << ", ";
    std::cout << ": \n";
    for (const auto &kt : siaVacsFull[is]) {
      auto iv = kt;
      if (iv < 0) continue;
      std::cout << "\t";
      for (const auto &jt : atoms[0][iv]) std::cout << jt << ", ";
      std::cout << "\n";
    }
    std::cout << "---\n";
  }
  std::cout << "vacSias before : \n";
  for (auto iv = 0; iv < vacSias.size(); iv++) {
    for (const auto &jt : atoms[0][iv]) std::cout << jt << ", ";
    std::cout << ": \n";
    for (const auto &kt : vacSias[iv]) {
      auto is = kt.second;
      if (is < 0) continue;
      std::cout << "\t";
      for (const auto &jt : atoms[1][is]) std::cout << jt << ", ";
      std::cout << "\n";
    }
    std::cout << "---\n";
  }
  */
  for (auto iv = 0; iv < vacSias.size(); iv++) {
    for (auto j = 1; j < vacSias[iv].size(); j++) {
      auto is = vacSias[iv][j].second;
      for (auto vac : siaVacsFull[is]) {
        if (vac == iv || !vacSias[vac].empty()) continue;
        auto disp = avi::calcDistSqr(atoms.sias[is].first, atoms.vacs[vac], box);
        vacSias[vac].push_back(std::make_pair(disp, is));
        vacSias[iv][j].second = -1;
      }
    }
  }
 
  /*
  std::cout << "atoms[0]: \n";
  for (auto it : atoms[0]) {
    for (auto jt : it) std::cout << jt << ", ";
    std::cout << '\n';
  }
  std::cout << "atoms[1]: \n";
  for (auto it : atoms[1]) {
    for (auto jt : it) std::cout << jt << ", ";
    std::cout << '\n';
  }
  auto temp = 0;
  std::cout << "vacSias: \n";
  for (auto it : vacSias) {
    for (auto jt : it) std::cout << jt.first << ", " << jt.second << "; ";
    std::cout << temp++ << "---\n";
  }
  std::cout << "vacSias after : \n";
  for (auto iv = 0; iv < vacSias.size(); iv++) {
    for (const auto &jt : atoms[0][iv]) std::cout << jt << ", ";
    std::cout << ": \n";
    for (const auto &kt : vacSias[iv]) {
      std::cout << "\t";
      auto is = kt.second;
      if (is < 0) continue;
      for (const auto &jt : atoms[1][is]) std::cout << jt << ", ";
      std::cout << "\n";
    }
    std::cout << "---\n";
  }
  */
  std::vector<int> latticeSiteGroups;
  //std::vector<size_t> ids; // TODO fill these
  std::vector<std::vector<double>> ids; // TODO fill these
  latticeSiteGroups.reserve(vacSias.size());
  for (auto i = 0; i < vacSias.size(); i++) {
    if (vacSias[i].size() == 1 && vacSias[i][0].first < threshDistSqr) continue;
    std::sort(begin(vacSias[i]), end(vacSias[i]));
    bool isAnnihilated = (!vacSias[i].empty() && vacSias[i][0].first < recombDistSqr);
    //std::cout << i << ": " << vacSias[i].size() << '\n' << std::flush;
    defects.emplace_back(atoms.vacs[i], false, defects.size() + 1, !isAnnihilated);  // vacancy
    if (!vacSias[i].empty()) {
      defects.emplace_back(atoms.sias[vacSias[i][0].second].first, true, defects.size() + 1, !isAnnihilated); // interstitial
      if (isIncludeId) ids.push_back(std::vector<double>{double(vacSias[i][0].second)});
      if (isExtraCol) {
        auto &dst = ids[ids.size()-1];
        auto &src = std::get<2>(statoms)[ids[ids.size()-1][0]];
        dst.insert(dst.end(), src.begin(), src.end());
      }
    }
    auto j = 1;
    for (; j < vacSias[i].size() && vacSias[i][j].first < vacGroupDistSqr; j++) {
      auto is = vacSias[i][j].second;
      if (is < 0) continue;
      defects.emplace_back(atoms.sias[is].first, true, defects.size() + 1, true); // interstitial
      if (isIncludeId) ids.push_back(std::vector<double>{double(is)});
      if (isExtraCol) {
        auto &dst = ids[ids.size()-1];
        auto &src = std::get<2>(statoms)[ids[ids.size()-1][0]];
        dst.insert(dst.end(), src.begin(), src.end());
      }
    }
    for (;j < vacSias[i].size(); j++) {
      auto is = vacSias[i][j].second;
      if (is < 0) continue;
      freeIs.push_back(is);
    }
    //std::cout << defects.size() << '\n';
    latticeSiteGroups.push_back(defects.size());
  }
  for (auto it : freeIs) {
    defects.emplace_back(atoms.sias[it].first, true, defects.size() + 1, true); // interstitial
    if (isIncludeId) ids.push_back(std::vector<double>{double(it)});
    if (isExtraCol) {
      auto &dst = ids[ids.size()-1];
      auto &src = std::get<2>(statoms)[ids[ids.size()-1][0]];
      dst.insert(dst.end(), src.begin(), src.end());
    }
  }
  /*
  std::cout << "latticeGroups: \n";
  for (auto it : latticeSiteGroups) {
    std::cout << it << ", ";
  }
  std::cout << '\n';
  */
  //std::cout <<"finnal defects size: " << defects.size() << '\n';
  return std::make_tuple(std::get<0>(statoms), avi::ErrorStatus::noError, std::move(defects), std::move(latticeSiteGroups), ids);
}
