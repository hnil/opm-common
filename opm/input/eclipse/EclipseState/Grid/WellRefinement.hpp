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
#ifndef OPM_WELL_REFINEMENT_HPP
#define OPM_WELL_REFINEMENT_HPP

#include <array>
#include <string>
#include <vector>

namespace Opm {

class Deck;
class DeckKeyword;

/// One WELLREF record: nested refinement rings around the wells matching a
/// pattern. Ring widths are coarse cells, outermost first; every ring is
/// refined NX x NY x NZ relative to the ring outside it.
struct WellRefinementRecord
{
    std::string pattern;
    std::array<int,3> factor{3, 3, 1};
    std::vector<int> rings;
    bool allLayers = false;     //!< every layer instead of the perforated ones
    int extraLayers = 0;        //!< layers above and below the perforated ones
};

/// A CARFIN block the expansion produces: 1-based inclusive cell range in the
/// parent grid's own cells, and the block's refined dimensions.
struct WellRefinementBox
{
    std::string name;
    std::string parent{"GLOBAL"};
    std::array<int,3> lo{};
    std::array<int,3> hi{};
    std::array<int,3> nxyz{};
};

std::vector<WellRefinementRecord> readWellRefinement(const DeckKeyword& keyword);

/// The CARFIN blocks for every WELLREF keyword in the deck, from the deck's
/// DIMENS, WELSPECS and COMPDAT records (wells with a trajectory keyword only
/// are skipped with a warning). Parents precede their children.
std::vector<WellRefinementBox> wellRefinementBoxes(const Deck& deck);

/// Ring construction on its own, for callers that have the seeds already:
/// one 0-based inclusive box per well, grid dimensions, the records.
std::vector<WellRefinementBox>
wellRefinementBoxes(const std::vector<WellRefinementRecord>& records,
                    const std::vector<std::pair<std::string, std::array<std::array<int,3>,2>>>& wellSeeds,
                    const std::array<int,3>& dims);

} // namespace Opm

#endif // OPM_WELL_REFINEMENT_HPP
