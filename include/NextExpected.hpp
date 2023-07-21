/*!
 * @file
 * helper functions for general usage
 * */

#ifndef NEXTEXPECTED_ANUVIKAR_HPP
#define NEXTEXPECTED_ANUVIKAR_HPP

#include <helper.hpp>

namespace avi {
struct NextExpected {
private:
  Coords _min, _max, _minCur, _maxCur, _maxFinal;
  Coords _cur;

public:
  NextExpected(Coords mn, Coords mx, Coords maxCur)
      : _min{mn}, _max{mx}, _cur{mn}, _minCur{mn}, _maxCur{maxCur}, _maxFinal{mx} {}

  NextExpected(Coords mn, Coords mx, Coords maxCur, Coords maxFinal)
      : _min{mn}, _max{mx}, _cur{mn}, _minCur{mn}, _maxCur{maxCur}, _maxFinal{maxFinal} {}

  NextExpected() = default;
  auto min() const { return _min; }
  auto max() const { return _max; }
  auto minCur() const { return _minCur; }
  auto maxCur() const { return _maxCur; }
  const Coords &cur() const { return _cur; }
  bool allMax() const {
    constexpr auto epsilon = 1e-4; // std::numeric_limits<double>::epsilon();
    for (size_t i = 0; i < _cur.size(); i++) {
      if ((_cur[i] + epsilon) < _maxFinal[i]) return false;
    }
    return true;
  }
  const Coords &increment();
  const Coords &incrementFcc();
  int mode{0};
};
} // namespace avi
#endif
