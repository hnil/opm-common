/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#define BOOST_TEST_MODULE WellRefinementTests

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/EclipseState/Grid/WellRefinement.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <array>
#include <string>

namespace {

const std::string deckText = R"(
RUNSPEC
DIMENS
 13 22 31 /
GRID
WELLREF
 'PROD1' 3 3 1 1 1 /
 'NOSUCH' 2 2 1 1 /
/
SCHEDULE
WELSPECS
 'PROD1' 'G' 6 3 1* 'OIL' /
 'INJ2'  'G' 12 20 1* 'WATER' /
/
COMPDAT
 'PROD1' 6 3 17 20 'OPEN' /
 'INJ2' 1* 1* 19 21 'OPEN' /
/
)";

std::array<int,3> ints(const Opm::DeckRecord& rec, const char* a, const char* b, const char* c)
{
    return { rec.getItem(a).get<int>(0), rec.getItem(b).get<int>(0), rec.getItem(c).get<int>(0) };
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(RingsBecomeNestedCarfinBlocks)
{
    const auto deck = Opm::Parser{}.parseString(deckText);

    // Two rings around PROD1 (6,3,17-20): the outer box, the inner nested in it.
    BOOST_REQUIRE_EQUAL(deck.count("CARFIN"), std::size_t{2});
    BOOST_CHECK_EQUAL(deck.count("ENDFIN"), std::size_t{2});

    const auto carfins = deck.getKeywordList("CARFIN");
    const auto& outer = carfins[0]->getRecord(0);
    BOOST_CHECK_EQUAL(outer.getItem("NAME").getTrimmedString(0), "WR1R1B1");
    BOOST_CHECK_EQUAL(outer.getItem("PARENT").getTrimmedString(0), "GLOBAL");
    BOOST_CHECK((ints(outer, "I1", "J1", "K1") == std::array<int,3>{4, 1, 16}));
    BOOST_CHECK((ints(outer, "I2", "J2", "K2") == std::array<int,3>{8, 5, 21}));
    BOOST_CHECK((ints(outer, "NX", "NY", "NZ") == std::array<int,3>{15, 15, 6}));

    // The inner ring is one coarse cell inside the outer in I, J and K, given in
    // the outer block's refined cells: coarse [5,7]x[2,4]x[17,20] ->
    // (1..3)*3 = 4..12 etc.
    const auto& inner = carfins[1]->getRecord(0);
    BOOST_CHECK_EQUAL(inner.getItem("NAME").getTrimmedString(0), "WR1R2B1");
    BOOST_CHECK_EQUAL(inner.getItem("PARENT").getTrimmedString(0), "WR1R1B1");
    BOOST_CHECK((ints(inner, "I1", "J1", "K1") == std::array<int,3>{4, 4, 2}));
    BOOST_CHECK((ints(inner, "I2", "J2", "K2") == std::array<int,3>{12, 12, 5}));
    BOOST_CHECK((ints(inner, "NX", "NY", "NZ") == std::array<int,3>{27, 27, 4}));

    // The blocks sit in the GRID section, right after WELLREF.
    std::size_t wellref = 0, first = 0, idx = 0;
    for (const auto& kw : deck) {
        if (kw.name() == "WELLREF") wellref = idx;
        if (kw.name() == "CARFIN" && first == 0) first = idx;
        ++idx;
    }
    BOOST_CHECK_EQUAL(first, wellref + 1);
}

BOOST_AUTO_TEST_CASE(DefaultedCompdatPositionTakesTheWellHead)
{
    const auto deck = Opm::Parser{}.parseString(deckText);
    const auto boxes = Opm::wellRefinementBoxes(deck);
    BOOST_REQUIRE_EQUAL(boxes.size(), std::size_t{2});

    // INJ2 is matched by no record, so nothing is generated for it, but its
    // defaulted COMPDAT I/J must resolve to the head without throwing.
    const Opm::WellRefinementRecord rec{ "INJ2", {2, 2, 1}, {1}, false, 0 };
    const auto one = Opm::wellRefinementBoxes(
        { rec }, { { "INJ2", { std::array<int,3>{11, 19, 18}, std::array<int,3>{11, 19, 20} } } },
        { 13, 22, 31 });
    BOOST_REQUIRE_EQUAL(one.size(), std::size_t{1});
    BOOST_CHECK((one[0].lo == std::array<int,3>{11, 19, 19}));
    BOOST_CHECK((one[0].hi == std::array<int,3>{13, 21, 21}));
    BOOST_CHECK((one[0].nxyz == std::array<int,3>{6, 6, 3}));
}

namespace {

const std::string nestedDeck = R"(RUNSPEC
DIMENS
5 5 1 /
GRID
CARFIN
'OUT'  2  4  2  4  1  1  9  9  1 /
ENDFIN
CARFIN
'IN'   4  6  4  6  1  1  9  9  1  1* 'OUT' /
ENDFIN
INIT
DX
 25*100 /
DY
 25*100 /
DZ
 25*10 /
TOPS
 25*2000 /
PORO
 25*0.3 /
PERMX
 25*500 /
PERMY
 25*500 /
PERMZ
 25*50 /
SCHEDULE
WELSPECL
 'PROD' 'G' 'IN' 5 5 1* 'OIL' /
/
COMPDATL
 'PROD' 'IN' 5 5 1 1 'OPEN' 1* 1* 0.2 /
/
)";

} // anonymous namespace

BOOST_AUTO_TEST_CASE(CompdatlInNestedLgr)
{
    const auto deck = Opm::Parser{}.parseString(nestedDeck);
    const auto es = Opm::EclipseState { deck };
    const auto labels = es.getInputGrid().get_all_labels();
    BOOST_REQUIRE_EQUAL(labels.size(), std::size_t{3});
    BOOST_CHECK_EQUAL(labels[1], "OUT");
    BOOST_CHECK_EQUAL(labels[2], "IN");

    const auto sched = Opm::Schedule { deck, es };
    const auto& well = sched.getWell("PROD", 0);
    BOOST_CHECK(well.is_lgr_well());
    BOOST_CHECK_EQUAL(well.get_lgr_well_tag().value(), "IN");
    const auto& conns = well.getConnections();
    BOOST_REQUIRE_EQUAL(conns.size(), std::size_t{1});
    BOOST_CHECK_EQUAL(conns[0].get_lgr_level(), 2);
    BOOST_CHECK_EQUAL(conns[0].getI(), 4);
    BOOST_CHECK_EQUAL(conns[0].getJ(), 4);
    BOOST_CHECK(conns[0].CF() > 0.0);
}
