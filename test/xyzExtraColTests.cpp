#include <catch.hpp>
#include <xyzReader.hpp>

using namespace avi;

TEST_CASE("Extra coordinates from a parcas xyz file line",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c{{0.,0.,0.}};
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordParcas("Fe  1 -76.770403   +7.2e2   .7",
                                               avi::frameStatus::inFrame, 3, 2, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == Approx(1.0));
    std::tie(ls, c, ec) =
        avi::getCoordParcas("what   +76.770403   -7.2e2   0.700 Frame",
                                 avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordParcas(
        "  what 2.000000 +3.2e2 4f 5.0 6.0 7.0 8.0 whatever", avi::frameStatus::inFrame, 2, 6, 7);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{2.0, 320, 4.0}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(6.0));
    CHECK(ec[1] == Approx(7.0));
    // garbage
    std::tie(ls, c, ec) =
        avi::getCoordParcas("r3   +76.770403   -7.2e2   0.700 Frame",
                                 avi::frameStatus::inFrame, 2, 5, 5);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) =
        avi::getCoordParcas("garbage 2.0 3.0 4.0 5.0 6.0 7.0", avi::frameStatus::inFrame, 2, 6, 8);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordParcas("-76.770403   +7.2e2   .7, 1.0",
                                               avi::frameStatus::inFrame, 1, 4, 4);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, +7.2e2, .7}});
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == Approx(1.0));
  }
}

TEST_CASE("Extra coordinates from a line from lammps xyz file",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c{{0.,0.,0.}};
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordLammps("Fe   -76.770403   +7.2e2   .7 5.0 6.0 7.0 8.0 what",
                                               avi::frameStatus::inFrame, 2, 6, 7);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(6.0));
    CHECK(ec[1] == Approx(7.0));
    std::tie(ls, c, ec) =
        avi::getCoordLammps("54   +76.770403   -7.2e2   0.700 ITEM:",
                                 avi::frameStatus::inFrame, 2, 1, 1);
    CHECK(ls == avi::lineStatus::coords);
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == Approx(54.0));
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "  what if 3 4 5 6 0.000000 +7.2e2 3f whatever  ", avi::frameStatus::inFrame, 7, 4, 5);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{0.0, 720, 3.0}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(4.0));
    CHECK(ec[1] == Approx(5.0));
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "no  what if 0.000000 +7.2e2 3f 2.0 ", avi::frameStatus::inFrame, 0, 4, 4);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{720, 3.0, 2.0}});
    REQUIRE(ec.size() == 0); // invalid precondition
    // garbage
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "no  what if 0.000000 +7.2e2 3f 2.0 whatever  ", avi::frameStatus::inFrame, 4, 7, 8);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) =
        avi::getCoordLammps("1.0 2.0 3.0 4.0", avi::frameStatus::inFrame, 1, 4, 6);
    CHECK(ls == avi::lineStatus::garbage);
  }
}

TEST_CASE("Extra coordinates from a line from Generic file",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c{{0.,0.,0.}};
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordGeneric("Fe   -76.770403   +7.2e2   .7 5.0 6.0 7.0 8.0 what",
                                               avi::frameStatus::inFrame, 2, 6, 7);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(6.0));
    CHECK(ec[1] == Approx(7.0));
    std::tie(ls, c, ec) =
        avi::getCoordGeneric("54   +76.770403   -7.2e2   0.700 ITEM:",
                                 avi::frameStatus::inFrame, 2, 1, 1);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == Approx(54.0));
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordGeneric(
        "  what if 3 4 5 6 0.000000 +7.2e2 3f whatever  ", avi::frameStatus::inFrame, 7, 4, 5);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    CHECK(c == Coords{{0.0, 720, 3.0}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(4.0));
    CHECK(ec[1] == Approx(5.0));
    std::tie(ls, c, ec) = avi::getCoordGeneric(
        "no  what if 0.000000 +7.2e2 3f 2.0 ", avi::frameStatus::inFrame, 0, 4, 4);
    CHECK(ls == avi::lineStatus::frameBorder);
    //CHECK(c == Coords{{0.0, 720, 3.0}});
    REQUIRE(ec.size() == 0); // invalid precondition
    // garbage
    std::tie(ls, c, ec) = avi::getCoordGeneric(
        "no  what if 0.000000 +7.2e2 3f 2.0 whatever  ", avi::frameStatus::inFrame, 4, 7, 8);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) =
        avi::getCoordGeneric("1.0 2.0 3.0 4.0", avi::frameStatus::inFrame, 1, 4, 6);
    CHECK(ls == avi::lineStatus::frameBorder);
  }
}

TEST_CASE("Extra coordinates from a line from Cdb file",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c{{0.,0.,0.}};
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordCdb("Fe   -76.770403   +7.2e2   .7 5.0 6.0 7.0 8.0 what",
                                               avi::frameStatus::inFrame, "Fe", 2, 6, 7);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == Approx(6.0));
    CHECK(ec[1] == Approx(7.0));
    std::tie(ls, c, ec) =
        avi::getCoordCdb("54   +76.770403   -7.2e2   0.700 ITEM:",
                                 avi::frameStatus::inFrame, "54", 2, 1, 1);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    REQUIRE(ec.size() == 0); // invalid case for cdb
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordCdb(
        "  what if 3 4 5 6 0.000000 +7.2e2 3f whatever  ", avi::frameStatus::inFrame, "what", 7, 4, 5);
    CHECK(ls == avi::lineStatus::inFrameCoords);
    CHECK(c == Coords{{0.0, 720, 3.0}});
    REQUIRE(ec.size() == 0); // invalid case for cdb
    std::tie(ls, c, ec) = avi::getCoordCdb(
        "no  what if 0.000000 +7.2e2 3f 2.0 ", avi::frameStatus::inFrame, "no", 0, 4, 4);
    CHECK(ls == avi::lineStatus::garbage);
    REQUIRE(ec.size() == 0);
    // garbage
    std::tie(ls, c, ec) = avi::getCoordCdb(
        "no  what if 0.000000 +7.2e2 3f 2.0 whatever  ", avi::frameStatus::inFrame, "no", 4, 7, 8);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) =
        avi::getCoordCdb("1.0 2.0 3.0 4.0", avi::frameStatus::inFrame, "1.0", 2, 5, 6);
    CHECK(ls == avi::lineStatus::frameBorder);
  }
          }

  TEST_CASE(
    "Read displaced (interstitial and vacancy pair) coordinates with extra cols from a line",
    "[defectsTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    std::array<Coords, 2> c;
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ");
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == std::array<Coords, 2>{{Coords{{-76.770403, 720, 0.7}},
                                      Coords{{76.770403, -720, 0.7}}}});
    CHECK(ec.empty());
    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, 1);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == std::array<Coords, 2>{{Coords{{-76.770403, 720, 0.7}},
                                      Coords{{76.770403, -720, 0.7}}}});
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == 6891.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == std::array<Coords, 2>{{Coords{{-76.770403, 720, 0.7}},
                                      Coords{{76.770403, -720, 0.7}}}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == 6891.0);
    CHECK(ec[1] == 112087.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 2, 3);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == std::array<Coords, 2>{{Coords{{-76.770403, 720, 0.7}},
                                      Coords{{76.770403, -720, 0.7}}}});
    REQUIRE(ec.size() == 2);
    CHECK(ec[0] == 112087.0);
    CHECK(ec[1] == 1.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, 3);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == std::array<Coords, 2>{{Coords{{-76.770403, 720, 0.7}},
                                      Coords{{76.770403, -720, 0.7}}}});
    REQUIRE(ec.size() == 3);
    CHECK(ec[0] == 6891.0);
    CHECK(ec[1] == 112087.0);
    CHECK(ec[2] == 1.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, 4);
    CHECK(ls == avi::lineStatus::garbage);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, -2);
    CHECK(ls == avi::lineStatus::coords);
    REQUIRE(ec.size() == 3);
    CHECK(ec[0] == 6891.0);
    CHECK(ec[1] == 112087.0);
    CHECK(ec[2] == 1.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", 1, -1);
    CHECK(ls == avi::lineStatus::coords);
    REQUIRE(ec.size() == 3);
    CHECK(ec[0] == 6891.0);
    CHECK(ec[1] == 112087.0);
    CHECK(ec[2] == 1.0);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 -76.770403   +7.2e2   .7  +76.770403   "
                                    "-7.2e2   0.700 6891 112087 1 ", -2, -2);
    CHECK(ls == avi::lineStatus::coords);
    REQUIRE(ec.size() == 3);
    CHECK(ec[0] == 6891.0);
    CHECK(ec[1] == 112087.0);
    CHECK(ec[2] == 1.0);

    std::tie(ls, c, ec) = avi::getCoordDisplaced(
        "1  -76.770403   +7.2e2  -76.770403   +7.2e2   .7  334 3 what", 1, 1);
    CHECK(ls == avi::lineStatus::coords);
    REQUIRE(ec.size() == 1);
    CHECK(ec[0] == 3.0);

    std::tie(ls, c, ec) = avi::getCoordDisplaced(
        "1  -76.770403   +7.2e2  -76.770403   +7.2e2   .7  334", 1, 1);
    CHECK(ls == avi::lineStatus::garbage);

    std::tie(ls, c, ec) =
        avi::getCoordDisplaced("1 34 2.5 3.5 1.0 3.0 3.5 what", 1, 1);
    CHECK(ls == avi::lineStatus::garbage);
  }
}