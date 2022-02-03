#include <catch.hpp>
#include <xyzReader.hpp>

using namespace avi;

TEST_CASE("Extra coordinates from a parcas xyz file line",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c;
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordParcas("Fe  1 -76.770403   +7.2e2   .7",
                                               avi::frameStatus::inFrame, 3, 2, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    std::tie(ls, c, ec) =
        avi::getCoordParcas("what   +76.770403   -7.2e2   0.700 Frame",
                                 avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordParcas(
        "  what 0.000000 +7.2e2 3f whatever  ", avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{0.0, 720, 3.0}});
    // garbage
    std::tie(ls, c, ec) =
        avi::getCoordParcas("53   +76.770403   -7.2e2   0.700 Frame",
                                 avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) =
        avi::getCoordParcas("garbage", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordParcas("-76.770403   +7.2e2   .7",
                                               avi::frameStatus::inFrame, 1);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, +7.2e2, .7}});
    std::tie(ls, c, ec) =
        avi::getCoordParcas("what 34 2.5", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordParcas(
        "-76.770403   +7.2e2   .7 garbage", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordParcas("what 34 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordParcas("what Frame 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    // border
    std::tie(ls, c, ec) = avi::getCoordParcas("Frame 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) =
        avi::getCoordParcas("Frame", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) = avi::getCoordParcas("  Frame 0.000000 +7.2e2 3f",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
  }
}

TEST_CASE("Extra coordinates from a line from lammps xyz file",
          "[extraColTest]") {
  SECTION("Normal cases") {
    avi::lineStatus ls;
    Coords c;
    std::vector<double> ec;
    // coords
    std::tie(ls, c, ec) = avi::getCoordLammps("Fe   -76.770403   +7.2e2   .7",
                                               avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{-76.770403, 720, 0.7}});
    std::tie(ls, c, ec) =
        avi::getCoordLammps("54   +76.770403   -7.2e2   0.700 ITEM:",
                                 avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{76.770403, -720, 0.7}});
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "  what if 0.000000 +7.2e2 3f whatever  ", avi::frameStatus::inFrame, 3);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{0.0, 720, 3.0}});
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "no  what if 0.000000 +7.2e2 3f 2.0 ", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::coords);
    CHECK(c == Coords{{720, 3.0, 2.0}});
    // garbage
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "no  what if 0.000000 +7.2e2 3f 2.0 whatever  ", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) =
        avi::getCoordLammps("garbage", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordLammps("-76.770403   +7.2e2   .7",
                                               avi::frameStatus::inFrame, 2);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) =
        avi::getCoordLammps("what 34 2.5", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordLammps(
        "-76.770403   +7.2e2   .7 garbage", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordLammps("what 34 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    std::tie(ls, c, ec) = avi::getCoordLammps("what ITEM: 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::garbage);
    // border
    std::tie(ls, c, ec) = avi::getCoordLammps("ITEM: 2.5 garbage",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) =
        avi::getCoordLammps("ITEM:", avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
    std::tie(ls, c, ec) = avi::getCoordLammps("  ITEM: 0.000000 +7.2e2 3f",
                                               avi::frameStatus::inFrame, 0);
    CHECK(ls == avi::lineStatus::frameBorder);
  }
}