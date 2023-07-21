#include <algorithm>

#include <AddOffset.hpp>
#include <optimizeOrigin.hpp>

// removes the overall offset not more than latConst / 4 from each coordinate of
// the origin
avi::Coords
avi::correctOrigin(const std::vector<avi::Coords> &atoms,
                        avi::Coords origin, double latConst) {
  std::vector<double> allOffsets;
  allOffsets.resize(atoms.size());
  const int topOffsetsToIgnore = atoms.size() > 2000 ? 1000 : atoms.size() / 10;
  avi::Coords res{{0.0, 0.0, 0.0}};
  for (int i = 0; i < 3; i++) {
    auto offsets = 0.0;
    for (size_t j = 0; j < atoms.size(); j++) {
      const auto &atom = atoms[j];
      auto offset = std::fmod(atom[i] - (origin[i] * latConst), latConst / 2.0);
      if (offset > latConst / 4.0) offset -= latConst / 2.0;
      allOffsets[j] = offset;
    }
    auto toConsider = allOffsets.size() - topOffsetsToIgnore;
    std::nth_element(
        begin(allOffsets), begin(allOffsets) + toConsider, end(allOffsets),
        [](double a, double b) { return std::fabs(a) < std::fabs(b); });
    for (size_t j = 0; j < toConsider; j++) {
      offsets += allOffsets[j];
    }
    res[i] = origin[i] + ((offsets / toConsider) / latConst);
  }
  return res;
}

avi::Coords
avi::estimateOrigin(const std::vector<avi::Coords> &atoms,
                         const double &latConst) {
  auto res = std::min_element(
      begin(atoms), end(atoms),
      [](const avi::Coords &a, const avi::Coords &b) {
        return (a[0] + a[1] + a[2]) < (b[0] + b[1] + b[2]);
      });
  auto originEstimated = *res;
  for (size_t i = 0; i < originEstimated.size(); i++)
    originEstimated[i] /= latConst;
  // we can return at this point but being extra careful with following
  // superflous steps
  auto obj2 = avi::AddOffset{latConst, "bcc", avi::StaticCoords{{originEstimated[0], originEstimated[1], originEstimated[2]}}};
  return correctOrigin(atoms, std::get<0>(obj2(*res)),
                       latConst); // rounded off and brought to first unit cell
}

/*
avi::Coords avi::estimateOrigin(const std::vector<avi::Coords>&
atoms, const avi::InputInfo& info) { auto res =
std::min_element(begin(atoms), end(atoms), [](const avi::Coords& a, const
avi::Coords& b) { return (a[0] + a[1] + a[2]) < (b[0] + b[1] + b[2]);
  });
  auto originEstimated = *res;
  for (size_t i = 0; i < originEstimated.size(); i++) originEstimated[i] /=
info.latticeConst;
  // we can return at this point but being extra careful with following
superflous steps auto obj2 = avi::AddOffset{info.latticeConst, "bcc",
originEstimated}; auto correct = correctOrigin(atoms,
std::get<0>(obj2(originEstimated)), info.latticeConst); // rounded off and
brought to first unit cell auto obj3 = avi::AddOffset{info.latticeConst,
"bcc", correct}; return std::get<0>(obj3(correct));
}
*/
