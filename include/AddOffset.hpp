/*!
 * @file
 * class for calculating offset & closest lattice site of an atom
 * */
#ifndef ADDOFFSET_ANUVIKAR_HPP
#define ADDOFFSET_ANUVIKAR_HPP

#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <helper.hpp>


namespace avi {

using offsetCoords = std::tuple<avi::Coords, double, avi::Coords, size_t>;

class AddOffset {
public:
  AddOffset(double latConst, std::string lattice, avi::StaticCoords origin);
  offsetCoords operator()(const avi::Coords &coords);

private:
  bool _isUnitcell(double x, double y, double z, double l,
                   std::array<double, 3> origin);
  void _bccUnitcell();
  void _fccUnitcell();
  double _calcDistMirror(std::array<long double, 3> a,
                         std::array<long double, 3> b, long double size,
                         std::array<int, 3> &mirror);

  std::vector<std::array<long double, 3>> _sites;
  std::array<long double, 3> _origin;
  long double _latConst;
  long double roundOffTo = 10000.0;
  size_t id {0};
};
} // namespace avi

#endif // ADDOFFSET_ANUVIKAR_HPP
