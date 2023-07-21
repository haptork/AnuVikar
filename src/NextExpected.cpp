#include <NextExpected.hpp>
#include <iostream>

const avi::Coords &avi::NextExpected::increment() {
  avi::Coords &c = _cur;
  constexpr auto epsilon = 1e-4;
  for (auto i : {2, 1, 0}) {
    //std::cout << "i is " << i << '\n';
    if (c[i] < (_maxCur[i] - epsilon)) { // normal coordinate increment
      if (i > 0) {
        c[i] += 1.0;
        //std::cout << "incrementing by one " << i << '\n';
      } else { // if (i == 0)
        c[i] += 0.5;
        //std::cout << "incrementing by half " << i << '\n';
        // ========== resetting min and max
        for (auto j = 0; j <= 2; ++j) {
          if (std::fabs(_minCur[j] - _min[j]) < epsilon) {
            _minCur[j] = _minCur[j] + 0.5;
          } else {
            _minCur[j] = _minCur[j] - 0.5;
          }
          if (std::fabs(_maxCur[j] - _max[j]) < epsilon) {
            _maxCur[j] = _maxCur[j] - 0.5;
          } else {
            _maxCur[j] = _maxCur[j] + 0.5;
          }
        }
        //std::cout << "new minmax\n";
        // ============
      }
      for (auto j = i + 1; j <= 2; ++j) {
        c[j] = _minCur[j];
      }
      return c;
    }
  }
  //std::cout << "outside\n";
  // TODO: enquire when does it reach here
  bool flag = false;
  for (int i = 0; i < _max.size(); i++) {
    if (_maxCur[i] != _max[i]) { flag = true; }
  }
  if (flag) {
    for (int i = 0; i < _max.size(); i++) {
      _maxCur[i] = _max[i];
      _minCur[i] = _min[i] + 0.5;
    }
    c[0] = _maxCur[0];
    c[1] = _minCur[1];
    c[2] = _minCur[2];
  }
  return c;
}

const avi::Coords &avi::NextExpected::incrementFcc() {
  avi::Coords &c = _cur;
  constexpr auto epsilon = 1e-4;
  //for (auto i : {2, 1, 0}) {
    if (c[2] < (_maxCur[2] - epsilon)) {
      c[2] += 1.0;
    } else if (c[1] < (_maxCur[1] - epsilon)) {
      c[1] += 0.5;
      int j = 2;
      if (std::fabs(_minCur[j] - _min[j]) < epsilon) {
        _minCur[j] = _minCur[j] + 0.5;
      } else {
        _minCur[j] = _minCur[j] - 0.5;
      }
      if (std::fabs(_maxCur[j] - _max[j]) < epsilon) {
        _maxCur[j] = _maxCur[j] - 0.5;
      } else {
        _maxCur[j] = _maxCur[j] + 0.5;
      }
      c[2] = _minCur[j];
    } else if (c[0] < (_maxCur[0] - epsilon)) {
      c[0] += 0.5;
      c[1] = _minCur[1];
      c[2] = _minCur[2];
    } else {
      //std::cout << "nothing";
  //std::cout << "outside\n";
  bool flag = false;
  for (int i = 0; i < _max.size(); i++) {
    if (_maxCur[i] != _max[i]) { flag = true; }
  }
  if (flag) {
    for (int i = 0; i < _max.size(); i++) {
      _maxCur[i] = _max[i];
      _minCur[i] = _min[i] + 0.5;
    }
    c[0] = _maxCur[0];
    c[1] = _minCur[1];
    c[2] = _minCur[2];
  }
  return c;

    }
    return c;
  //}
  // TODO: enquire when does it reach here
}
