#include <catch.hpp>

#include <AddOffset.hpp>
#include <NextExpected.hpp>
#include <results.hpp>
#include <xyz2defects.hpp>
#include <xyzReader.hpp>

#include <iostream>

using namespace avi;
SCENARIO("FCC- Find nearest lattice site for a coordinate given lattice structure - "
         "Addoffset",
         "[defectsTestFcc]") {
  SECTION("Normal Cases - With 0.0 as origin") {
    auto origin = Coords{{0.0, 0.0, 0.0}};
    // Case 1
    AddOffset a{1.0, "fcc", origin};
    auto res = a(origin);
    REQUIRE(std::get<0>(res)[0] == Approx(0.0));
    REQUIRE(std::get<0>(res)[1] == Approx(0.0));
    REQUIRE(std::get<0>(res)[2] == Approx(0.0));
    REQUIRE(std::get<1>(res) == Approx(0.0));
    REQUIRE(std::get<2>(res) == origin);
    // Case a
    auto resa = a(Coords{{0.4, 0.5, 0.5}});
    REQUIRE(std::get<0>(resa)[0] == Approx(0.0));
    REQUIRE(std::get<0>(resa)[1] == Approx(0.5));
    REQUIRE(std::get<0>(resa)[2] == Approx(0.5));
    REQUIRE(std::get<1>(resa) == Approx(0.4));
    // Case b
    auto cb = Coords{{2.0, 1.5, 3.71}};
    auto resb = a(cb);
    CHECK(std::get<0>(resb)[0] == Approx(2.0));
    CHECK(std::get<0>(resb)[1] == Approx(1.5));
    CHECK(std::get<0>(resb)[2] == Approx(3.5));
    REQUIRE(std::get<1>(resb) == Approx(0.21));
    // checking if the lattice site found is nearer than other sites around
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 1.5, 3.5}})); // 0.583
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 2.0, 4.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.5, 1.5, 4.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 1.0, 4.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{1.5, 1.5, 4.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.5, 1.5, 3.5}})); // 0.735
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{1.5, 2.0, 3.5}})); // 0.735
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{1.5, 1.0, 3.5}})); // 0.735
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.5, 2.0, 3.5}})); // 0.735
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 1.5, 4.0}})); // 0.800
    // Case c
    auto cc = Coords{{-0.1, -0.3, -0.71}};
    auto resc = a(cc);
    CHECK(std::get<0>(resc)[0] == Approx(-0.0));
    CHECK(std::get<0>(resc)[1] == Approx(-0.5));
    CHECK(std::get<0>(resc)[2] == Approx(-0.5));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{0.0, 0.0, 0.0}}));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{-0.5, -0.5, 0.0}}));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{-0.5, -0.5, -1.0}}));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{0.0, 0.0, -1.0}}));
    CHECK(std::get<1>(resc) == Approx(0.307).epsilon(0.01));
    // Case d
    auto cd = Coords{{-0.3, -0.3, -0.71}};
    auto resd = a(cd);
    CHECK(std::get<0>(resd)[0] == Approx(-0.5));
    CHECK(std::get<0>(resd)[1] == Approx(-0.5));
    CHECK(std::get<0>(resd)[2] == Approx(-1.0));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{0.0, 0.0, 0.0}}));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{-1.0, -1.0, -1.0}}));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{0.0, 0.0, -1.0}}));
    CHECK(std::get<1>(resd) >= calcDist(cd, Coords{{-0.5, -0.5, -0.5}}));
    CHECK(std::get<1>(resd) == Approx(0.405).epsilon(0.001));
  }
  SECTION("Normal Cases - With 0.25 as origin") {
    auto origin = Coords{{0.25, 0.25, 0.25}};
    // Case 1
    AddOffset a{1.0, "fcc", origin};
    auto res = a(origin);
    REQUIRE(std::get<0>(res)[0] == Approx(origin[0]));
    REQUIRE(std::get<0>(res)[1] == Approx(origin[1]));
    REQUIRE(std::get<0>(res)[2] == Approx(origin[2]));
    REQUIRE(std::get<1>(res) == Approx(0.0));
    REQUIRE(std::get<2>(res) == origin);
    // Case a
    auto resa = a(Coords{{0.64, 0.75, 0.75}});
    REQUIRE(std::get<0>(resa)[0] == Approx(0.25));
    REQUIRE(std::get<0>(resa)[1] == Approx(0.75));
    REQUIRE(std::get<0>(resa)[2] == Approx(0.75));
    REQUIRE(std::get<1>(resa) == Approx(0.39));
    // Case b
    auto cb = Coords{{2.0, 1.5, 3.71}};
    auto resb = a(cb);
    CHECK(std::get<0>(resb)[0] == Approx(1.75));
    CHECK(std::get<0>(resb)[1] == Approx(1.25));
    CHECK(std::get<0>(resb)[2] == Approx(3.75));
    // checking if the lattice site found is nearer than other sites around
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.25, 1.25, 3.25}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.25, 2.25, 3.25}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.25, 2.25, 4.25}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.25, 1.25, 4.25}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.75, 1.75, 3.75}}));
    // Case c
    auto cc = Coords{{-0.1, -0.3, -0.71}};
    auto resc = a(cc);
    CHECK(std::get<0>(resc)[0] == Approx(-0.25));
    CHECK(std::get<0>(resc)[1] == Approx(-0.25));
    CHECK(std::get<0>(resc)[2] == Approx(-0.75));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{0.25, 0.25, 0.25}}));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{-0.75, -0.75, -0.75}}));
    CHECK(std::get<1>(resc) <= calcDist(cc, Coords{{-0.25, -0.25, -1.25}}));
    CHECK(std::get<1>(resc) == Approx(0.163).epsilon(0.01));
    // Case d
    auto cd = Coords{{-0.3, -0.3, -0.71}};
    auto resd = a(cd);
    CHECK(std::get<0>(resd)[0] == Approx(-0.25));
    CHECK(std::get<0>(resd)[1] == Approx(-0.25));
    CHECK(std::get<0>(resd)[2] == Approx(-0.75));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{0.25, 0.25, 0.25}}));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{-0.25, -0.25, -0.25}}));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{-0.75, -0.75, -0.75}}));
    CHECK(std::get<1>(resd) <= calcDist(cd, Coords{{-0.25, -0.25, -1.25}}));
    CHECK(std::get<1>(resd) == Approx(0.081).epsilon(0.01));
  }
  SECTION("Normal Cases - With 0.5 as origin and Fe lattice constant (2.85)") {
    // TODO
    /*
    auto origin = Coords{{0.5, 0.5, 0.5}};
    auto originMult = Coords{{0.5 * 2.85, 0.5 * 2.85, 0.5 * 2.85}};
    // Case 1
    AddOffset a{2.85, "bcc", origin};
    auto res = a(originMult);
    CHECK(std::get<0>(res)[0] == Approx(origin[0]));
    CHECK(std::get<0>(res)[1] == Approx(origin[1]));
    CHECK(std::get<0>(res)[2] == Approx(origin[2]));
    CHECK(std::get<1>(res) == Approx(0.0));
    CHECK(std::get<2>(res) == originMult);
    // Case a
    auto resa = a(Coords{{2.83, 2.85, 2.85}});
    REQUIRE(std::get<0>(resa)[0] == Approx(1.0));
    REQUIRE(std::get<0>(resa)[1] == Approx(1.0));
    REQUIRE(std::get<0>(resa)[2] == Approx(1.0));
    REQUIRE(std::get<1>(resa) == Approx(0.02));
    // Case b
    auto cb = Coords{{5.79, 3.98, 7.71}};
    auto resb = a(cb);
    CHECK(std::get<0>(resb)[0] == Approx(2.00));
    CHECK(std::get<0>(resb)[1] == Approx(1.00));
    CHECK(std::get<0>(resb)[2] == Approx(3.00));
    // checking if the lattice site found is nearer than other sites around
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{1.5, 1.5, 2.5}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 2.00, 3.00}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 2.0, 3.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.0, 1.0, 4.0}}));
    CHECK(std::get<1>(resb) <= calcDist(cb, Coords{{2.5, 1.5, 3.5}}));
    */
  }
  SECTION("Edge Cases") {
    auto origin = Coords{{0.0, 0.0, 0.0}};
    // Case 1
    AddOffset a{1.0, "bcc", origin};
    auto c1 = Coords{{-0.0, -0.00001, 0.00001}};
    auto res = a(c1);
    REQUIRE(std::get<0>(res)[0] == Approx(0.0));
    REQUIRE(std::get<0>(res)[1] == Approx(0.0));
    REQUIRE(std::get<0>(res)[2] == Approx(0.0));
    REQUIRE(std::get<2>(res) == c1);
  }
}

SCENARIO("Fcc- Enumerate all lattice sites for a bcc in ascending order given min "
         "and max lattice site. The lattice site are relative sites (i.e. "
         "latticeConst does not affect the values.) The min and max given "
         "should also be relative. The simulation codes like lammps / parcas "
         "build bcc lattice such that both the intertwined simple cubic "
         "lattices in bcc have same number of unit-cells. The code can assume "
         "that min and max are given accordingly in valid forms. The code "
         "should accomodate for change in origin according to the min value. - "
         "NextExptected",
         "[defectsTestFcc]") {
  SECTION("Normal Case - 1") {
    // case 1
    auto origin = avi::Coords{{0., 0., 0.}};
    auto max = avi::Coords{{2.5, 2.5, 2.5}};
    auto maxInitial = getInitialMaxFcc(origin, max);
    NextExpected ne{
        origin, max,
        maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.cur() == Coords{{0, 0, 0}});
    CHECK(ne.incrementFcc() == Coords{{0, 0, 1}});
    CHECK(ne.incrementFcc() == Coords{{0, 0, 2}});
    CHECK(ne.incrementFcc() == Coords{{0, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 0.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 0.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 1, 0}});
    CHECK(ne.incrementFcc() == Coords{{0, 1, 1}});
    CHECK(ne.incrementFcc() == Coords{{0, 1, 2}});
    CHECK(ne.incrementFcc() == Coords{{0, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 1.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 1.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 2, 0}});
    CHECK(ne.incrementFcc() == Coords{{0, 2, 1}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0, 2, 2}});
    CHECK(ne.incrementFcc() == Coords{{0, 2.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 2.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0, 2.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 2, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 2, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 2, 2.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 2.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 2.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 2.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{1, 0, 0}});
    CHECK(ne.incrementFcc() == Coords{{1, 0, 1}});
    CHECK(ne.incrementFcc() == Coords{{1, 0, 2}});
    CHECK(ne.incrementFcc() == Coords{{1, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 0.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 0.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 1, 0}});
    CHECK(ne.incrementFcc() == Coords{{1, 1, 1}});
    CHECK(ne.incrementFcc() == Coords{{1, 1, 2}});
    CHECK(ne.incrementFcc() == Coords{{1, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 1.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 1.5, 2.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1, 2, 0}});
    CHECK(ne.incrementFcc() == Coords{{1, 2, 1}});
    CHECK(ne.incrementFcc() == Coords{{1, 2, 2}});
    CHECK(ne.incrementFcc() == Coords{{1, 2.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 2.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1, 2.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 2, 0.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.5, 2, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 2, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 2.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 2.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 2.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{2, 0, 0}});
    CHECK(ne.incrementFcc() == Coords{{2, 0, 1}});
    CHECK(ne.incrementFcc() == Coords{{2, 0, 2}});
    CHECK(ne.incrementFcc() == Coords{{2, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 0.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 0.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 1, 0}});
    CHECK(ne.incrementFcc() == Coords{{2, 1, 1}});
    CHECK(ne.incrementFcc() == Coords{{2, 1, 2}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 1.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 1.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 2, 0}});
    CHECK(ne.incrementFcc() == Coords{{2, 2, 1}});
    CHECK(ne.incrementFcc() == Coords{{2, 2, 2}});
    CHECK(ne.incrementFcc() == Coords{{2, 2.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 2.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{2, 2.5, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 0.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 1.5, 2}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 2, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 2, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.5, 2, 2.5}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 2.5, 0}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 2.5, 1}});
    CHECK(ne.incrementFcc() == Coords{{2.5, 2.5, 2}});
    REQUIRE(ne.allMax());
  }
  SECTION("Normal Case - 2") {
    // case 2
    auto origin = avi::Coords{{0.25, 0.25, 0.25}};
    auto max = avi::Coords{{2.75, 2.75, 2.75}};
    auto maxInitial = getInitialMaxFcc(origin, max);
    NextExpected ne{
        origin, max,
        maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.25
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.cur() == Coords{{0.25, 0.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 0.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 0.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 0.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 0.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 0.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 1.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 2.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 2.25, 1.25}});
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.incrementFcc() == Coords{{0.25, 2.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{0.25, 2.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 2.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.25, 2.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 0.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 1.75, 2.25}});
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.incrementFcc() == Coords{{0.75, 2.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 2.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 2.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{0.75, 2.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 2.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{0.75, 2.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 0.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 1.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 2.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 2.25, 1.25}});
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.incrementFcc() == Coords{{1.25, 2.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.25, 2.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 2.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.25, 2.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 0.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 1.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{1.75, 2.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 0.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 0.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 0.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 0.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 0.75, 1.75}});
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.incrementFcc() == Coords{{2.25, 0.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 1.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.25, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.25, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.25, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.75, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.75, 1.75}});
CHECK(ne.incrementFcc() == Coords{{2.25, 2.75, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 0.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.25, 1.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 1.75, 2.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 2.25, 0.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 2.25, 1.75}});
    REQUIRE_FALSE(ne.allMax());
CHECK(ne.incrementFcc() == Coords{{2.75, 2.25, 2.75}});
CHECK(ne.incrementFcc() == Coords{{2.75, 2.75, 0.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 2.75, 1.25}});
CHECK(ne.incrementFcc() == Coords{{2.75, 2.75, 2.25}});
    REQUIRE(ne.allMax());
  }
  SECTION("Edge Cases - 1") {
    // case 2
    /*
    auto origin = avi::Coords{{0.5, 0.5, 0.5}};
    auto max = avi::Coords{{2.00, 2.00, 2.00}};
    auto maxInitial = getInitialMax(origin, max);
    NextExpected ne{origin, max, maxInitial};
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.0, 1.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 1.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 2.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 2.0, 2.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.0, 1.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 1.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 2.0, 1.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.0, 2.0, 2.0}});
    REQUIRE(ne.allMax());
    */
  }
  SECTION("Edge Cases - Almost invalid input - 1") {
    // case 1
    /*
    auto origin = avi::Coords{{0.0, 0.0, 0.0}};
    auto max = avi::Coords{
        {2.00, 2.00, 2.00}}; // a valid max should have been 2.5, orign 0.0
    auto maxInitial = getInitialMax(origin, max);
    NextExpected ne{origin, max, maxInitial};
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.cur() == Coords{{0.0, 0.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 0.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 0.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 1.0, 0.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.0, 1.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 1.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 2.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 2.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{0.0, 2.0, 2.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{0.5, 0.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 0.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{0.5, 1.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.0, 0.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 0.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 0.0, 2.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.0, 1.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 1.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 1.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 2.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 2.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{1.0, 2.0, 2.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 0.5, 1.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{1.5, 1.5, 1.5}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.0, 0.5, 0.5}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 0.5, 1.5}});
    */
    /*
    CHECK(ne.incrementFcc() == Coords{{2.0, 0.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 0.0, 2.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 1.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 1.0, 1.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 1.0, 2.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.0, 2.0, 0.0}});
    CHECK(ne.incrementFcc() == Coords{{2.0, 2.0, 1.0}});
    REQUIRE_FALSE(ne.allMax());
    CHECK(ne.incrementFcc() == Coords{{2.0, 2.0, 2.0}});
    //REQUIRE(ne.max() ==  Coords{{0.0,0.0,0.0}});
    REQUIRE(ne.allMax());
    */
  }
}

 SCENARIO("Fcc- Given xyz coordinates of all the lattice atoms, output only the "
          "defects, labelled by interstitial and vacancy, annihilated (psuedo) "
          "or actual (true)",
          "[defectsTestFcc]") {
   SECTION("Normal cases") {
     SECTION("Single Dumbbell") {
       std::vector<std::tuple<Coords, double, Coords>> atoms;
       auto origin = avi::Coords{{0.25, 0.25, 0.25}};
       //auto max = avi::Coords{{10.75, 10.75, 10.75}};
       auto max = avi::Coords{{3.75, 3.75, 3.75}};
       auto maxInitial = getInitialMaxFcc(origin, max);
       NextExpected ne{
           origin, max,
           maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
       auto latticeConst = 1.00;
       AddOffset addOffset{latticeConst, "fcc", origin}; // Fe
       Coords interstitialCoord, vacancyCoord;
       auto i = 0;
       auto pickAt = 30; // this is ~ 10th atom
       //std::cout << "begins\n";
       while (true) {
         Coords c = ne.cur();
         // if (i != pickAt) {
         for (auto &jt : c) {
           jt *= latticeConst;
           jt += (i++ & 1) ? -0.05 : 0.05;
         }
         //}
         if (i == pickAt) {
           vacancyCoord = ne.cur();
           for (auto &jt : vacancyCoord)
             jt *= latticeConst;
           auto picked = atoms[pickAt / 5];
           for (int k = 0; k < 3; ++k) {
             auto diff = std::get<0>(picked)[k] - std::get<2>(picked)[k];
             c[k] = std::get<0>(picked)[k] - diff + 0.001;
           }
           interstitialCoord = c;
         }
         //std::cout << c[0] << ", " << c[1] << ", " << c[2] << '\n';
         atoms.emplace_back(addOffset(c));
         if (ne.allMax()) break;
         ne.incrementFcc();
       }
       InputInfo info;
       Config config;
       config.safeRunChecks = false;
       config.isIgnoreBoundaryDefects = false;
       info.latticeConst = latticeConst;
       info.originX = origin[0];
       info.originY = origin[1];
       info.originZ = origin[2];
       ExtraInfo extraInfo;
       auto fsAtoms = std::make_pair(avi::xyzFileStatus::reading, atoms);
       auto ungroupedDefectsDumbbellPair = atoms2defectsFcc(fsAtoms, info, extraInfo, config);
       auto ungroupedDefects = std::get<2>(ungroupedDefectsDumbbellPair);
       REQUIRE(ungroupedDefects.size() == 4); // 2 interstitials, 2 vacancies
       SECTION("Check cluster grouping") {
         int nDefects;
         double inClusterFractionI, inClusterFractionV;
         auto defects = groupDefects(ungroupedDefects, latticeConst);
         auto clusterSizeMap = clusterSizes(defects);
         REQUIRE(clusterSizeMap.size() ==
                 2); // one cluster of three and other of one
         SECTION("Check ndefects and cluster sizes") {
           std::tie(nDefects, inClusterFractionI, inClusterFractionV) =
               avi::getNDefectsAndClusterFractions(defects);
           REQUIRE(nDefects == 1);
           REQUIRE(inClusterFractionI == Approx(100.0));
           REQUIRE(inClusterFractionV == Approx(100.0));
           ignoreSmallClusters(defects, clusterSizeMap);
           std::tie(nDefects, inClusterFractionI, inClusterFractionV) =
               avi::getNDefectsAndClusterFractions(defects);
           REQUIRE(nDefects == 1);
           REQUIRE(inClusterFractionI == Approx(0.0));
           REQUIRE(inClusterFractionV == Approx(0.0)); // changed
           auto clusterIdMap = avi::clusterMapping(defects);
           REQUIRE(clusterIdMap.size() == 0); // 1 dumbbell
           auto clusterIVMap =
               avi::clusterIVType(clusterIdMap, clusterSizeMap);
           REQUIRE(clusterIVMap.size() == 0);
           int maxClusterSizeI, maxClusterSizeV;
           std::tie(maxClusterSizeI, maxClusterSizeV) =
               avi::getMaxClusterSizes(clusterSizeMap, clusterIdMap);
           REQUIRE(maxClusterSizeI == 0);
           REQUIRE(maxClusterSizeV == 0);
           /*
         for (auto x : defects) {
           for (auto c : std::get<0>(x)) std::cout << c << ", ";
           std::cout << std::get<1>(x) << ", " << std::get<2>(x) << ", " << std::get<3>(x);
           std::cout << std::endl;
           }
         */
           SECTION("Check cluster features") {
             auto feats = avi::clusterFeatures(
                 defects, clusterIdMap, clusterSizeMap, latticeConst);
             REQUIRE(feats.size() == 0);
             /*
             const auto &distFeat = std::get<0>(std::begin(feats)->second);
             REQUIRE(distFeat[0] + distFeat[distFeat.size() - 1] == 1.0); // TODO
             REQUIRE(distFeat[0] == Approx(2.0 / 6.0)); // TODO
             // REQUIRE(distFeat[distFeat.size() - 1] == Approx(4.0 / 6.0));
             REQUIRE(std::accumulate(begin(distFeat), end(distFeat), 0.0) ==
                     Approx(1.0));
             const auto &angleFeat = std::get<1>(std::begin(feats)->second);
             // REQUIRE(angleFeat[0] + angleFeat[angleFeat.size() - 1] == Approx(1.0));
             REQUIRE(std::accumulate(begin(angleFeat), end(angleFeat), 0.0) ==
                     Approx(1.0));
             const auto &adjFeat = std::get<2>(std::begin(feats)->second);
             REQUIRE(adjFeat[2] == Approx(1.0));
             REQUIRE(std::accumulate(begin(adjFeat), end(adjFeat), 0.0) ==
                     Approx(1.0));
                     */
           }
         }   // ndefects and cluster sizes
       }     // cluster grouping
     }       // End of Single Dumbbell test
    SECTION("big interstitial cluster") {
      std::vector<std::tuple<Coords, double, Coords>> atoms;
      auto origin = avi::Coords{{0.5, 0.5, 0.5}};
      auto max = avi::Coords{{6.0, 6.0, 6.0}};
      auto maxInitial = getInitialMaxFcc(origin, max);
      NextExpected ne{
          origin, max,
          maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
      auto latticeConst = 3.165; // W
      AddOffset addOffset{latticeConst, "fcc", origin};
      Coords lastInterstitialCoord, lastVacancyCoord;
      auto i = 0;
      auto pickAt = 300; // this is ~ 100th atom
      while (true) {
        // std::cout << ne.cur()[0] << ", " << ne.cur()[1] << ", " <<
        // ne.cur()[2] << " | " << ne.minCur()[0] << ", " << ne.maxCur()[0] <<
        // '\n';
        Coords c = ne.cur();
        // if (i != pickAt) {
        for (auto &jt : c) {
          jt *= latticeConst;
          jt += (i++ & 1) ? -0.11 : 0.11;
        }
        //}
        if (i % pickAt == 0) {
          lastVacancyCoord = ne.cur();
          for (auto &jt : lastVacancyCoord)
            jt *= latticeConst;
          auto picked = atoms[pickAt / 5];
          for (int k = 0; k < 3; ++k) {
            auto diff = std::get<0>(picked)[k] - std::get<2>(picked)[k];
            c[k] = std::get<0>(picked)[k] - diff + (0.001 * (i / pickAt));
          }
          lastInterstitialCoord = c;
        }
        atoms.emplace_back(addOffset(c));
        if (ne.allMax()) break;
        ne.incrementFcc();
      }
      InputInfo info;
      Config config;
      config.safeRunChecks = false;
      config.isIgnoreBoundaryDefects = false;
      info.latticeConst = latticeConst;
      info.originX = origin[0];
      info.originY = origin[1];
      info.originZ = origin[2];
      ExtraInfo extraInfo;
      auto fsAtoms = std::make_pair(avi::xyzFileStatus::reading, atoms);
      auto ungroupedDefectsDumbbellPair = atoms2defectsFcc(fsAtoms, info, extraInfo, config);
      auto ungroupedDefects = std::get<2>(ungroupedDefectsDumbbellPair);
      REQUIRE(ungroupedDefects.size() == 16);
      int nDefects;
      double inClusterFractionI, inClusterFractionV;
      std::tie(nDefects, inClusterFractionI, inClusterFractionV) =
          avi::getNDefectsAndClusterFractions(ungroupedDefects);
      SECTION("Check cluster grouping") {
        auto defects = groupDefects(ungroupedDefects, latticeConst);
        auto clusterSizeMap = clusterSizes(defects);
        REQUIRE(clusterSizeMap.size() == 7);
        SECTION("Check ndefects and cluster sizes") {
          std::tie(nDefects, inClusterFractionI, inClusterFractionV) =
              avi::getNDefectsAndClusterFractions(defects);
          REQUIRE(nDefects == 8);
          REQUIRE(inClusterFractionI == Approx(100.0));
          REQUIRE(inClusterFractionV == Approx(75.0));
          ignoreSmallClusters(defects, clusterSizeMap);
          std::tie(nDefects, inClusterFractionI, inClusterFractionV) =
              avi::getNDefectsAndClusterFractions(defects);
          REQUIRE(nDefects == 8);
          REQUIRE(inClusterFractionI == Approx(100.0));
          REQUIRE(inClusterFractionV == Approx(0.0)); // changed
          auto clusterIdMap = avi::clusterMapping(defects);
          REQUIRE(clusterIdMap.size() == 1); // 1 interstitial cluster
          REQUIRE(std::begin(clusterIdMap)->second.size() == 10);
          auto clusterIVMap =
              avi::clusterIVType(clusterIdMap, clusterSizeMap);
          REQUIRE(clusterIVMap.size() == 1);
          REQUIRE(std::begin(clusterIVMap)->second == 8); // surviving
          int maxClusterSizeI, maxClusterSizeV;
          std::tie(maxClusterSizeI, maxClusterSizeV) =
              avi::getMaxClusterSizes(clusterSizeMap, clusterIdMap);
          REQUIRE(maxClusterSizeI == 8);
          REQUIRE(maxClusterSizeV == 0);
        } // ndefects and cluster sizes
      }   // cluster grouping
    }
  }
  SECTION("Invalid input cases") {
    SECTION(
        "One atom kicked out of the box : unwrapped coordinates are invalid") {
      std::vector<std::tuple<Coords, double, Coords>> atoms;

      auto origin = avi::Coords{{0.25, 0.25, 0.25}};
      auto max = avi::Coords{{10.75, 10.75, 10.75}};
      auto maxInitial = getInitialMaxFcc(origin, max);
      NextExpected ne{
          origin, max,
          maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
      REQUIRE_FALSE(ne.allMax());
      auto latticeConst = 2.85;
      AddOffset addOffset{latticeConst, "fcc", origin}; // Fe
      auto i = 0;
      auto pickAt = 30; // this is ~ 10th atom
      while (true) {
        Coords c = ne.cur();
        // if (i != pickAt) {
        for (auto &jt : c) {
          jt *= latticeConst;
          jt += (i++ & 1) ? -0.15 : 0.15;
        }
        //}
        if (i == pickAt) {
          for (int k = 0; k < 3; ++k) {
            c[k] += latticeConst * 11;
          }
        }
        atoms.emplace_back(addOffset(c));
        if (ne.allMax()) break;
        ne.incrementFcc();
      }
      InputInfo info;
      Config config;
      config.safeRunChecks = false;
      info.latticeConst = latticeConst;
      info.originX = origin[0];
      info.originY = origin[1];
      info.originZ = origin[2];
      ExtraInfo extraInfo;
      auto fsAtoms = std::make_pair(avi::xyzFileStatus::reading, atoms);
      auto ungroupedDefectsDumbbellPair = atoms2defectsFcc(fsAtoms, info, extraInfo, config);
      auto ungroupedDefects = std::get<2>(ungroupedDefectsDumbbellPair);
      // it should have been 4 but now alot more defects are counted as the
      // box size is inferred from the atom coordinates, including the atom
      // that we kicked out of the box
      REQUIRE(ungroupedDefects.size() > 4);
      REQUIRE(ungroupedDefects.size() >
              400); // even greater than 100 times of actual
      auto interstitialCount = 0;
      auto vacancyCount = 0;
      for (const auto &it : ungroupedDefects) {
        if (DefectTWrap::isInterstitial(it))
          interstitialCount++;
        else
          vacancyCount++;
      }
      REQUIRE(vacancyCount > interstitialCount);
      //REQUIRE(interstitialCount == 0); 
      // TODO
    }
  }
  SECTION("Edge cases") {
    SECTION("Perfect lattice") {
      std::vector<std::tuple<Coords, double, Coords>> atoms;
      auto origin = avi::Coords{{0.0, 0.0, 0.0}};
      auto max = avi::Coords{{1.5, 1.5, 1.5}};
      auto maxInitial = getInitialMaxFcc(origin, max);
      NextExpected ne{
          origin, max,
          maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
      AddOffset addOffset{1.0, "fcc", origin};
      while (true) {
        atoms.emplace_back(addOffset(ne.cur()));
        if (ne.allMax()) break;
        ne.incrementFcc();
      }
      InputInfo info;
      Config config;
      config.safeRunChecks = false;
      info.latticeConst = 1.0;
      info.originX = origin[0];
      info.originY = origin[1];
      info.originZ = origin[2];
      ExtraInfo extraInfo;
      auto fsAtoms = std::make_pair(avi::xyzFileStatus::reading, atoms);
      auto ungroupedDefectsDumbbellPair = atoms2defectsFcc(fsAtoms, info, extraInfo, config);
      REQUIRE(std::get<2>(ungroupedDefectsDumbbellPair).empty());
    }
    SECTION("slightly shaken lattice") {
      std::vector<std::tuple<Coords, double, Coords>> atoms;
      auto origin = avi::Coords{{0.0, 0.0, 0.0}};
      auto max = avi::Coords{{1.5, 1.5, 1.5}};
      auto maxInitial = getInitialMaxFcc(origin, max);
      NextExpected ne{
          origin, max,
          maxInitial, maxInitial}; // The min and max for two unit-cells with origin 0.0
      AddOffset addOffset{1.0, "fcc", Coords{0.0, 0.0, 0.0}};
      auto i = 0;
      while (true) {
        Coords c = ne.cur();
        for (auto &jt : c) {
          jt += (i++ & 1) ? -0.1 : 0.15;
        }
        atoms.emplace_back(addOffset(c));
        if (ne.allMax()) break;
        ne.incrementFcc();
      }
      InputInfo info;
      ExtraInfo extraInfo;
      Config config;
      config.safeRunChecks = false;
      info.latticeConst = 1.0;
      info.originX = origin[0];
      info.originY = origin[1];
      info.originZ = origin[2];
      auto fsAtoms = std::make_pair(avi::xyzFileStatus::reading, atoms);
      auto ungroupedDefectsDumbbellPair = atoms2defectsFcc(fsAtoms, info, extraInfo, config);
      REQUIRE(std::get<2>(ungroupedDefectsDumbbellPair).empty());
    }
  }
}
