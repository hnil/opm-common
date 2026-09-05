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
#include <opm/input/eclipse/EclipseState/Grid/WellRefinement.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/shmatch.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace Opm {

std::vector<WellRefinementRecord> readWellRefinement(const DeckKeyword& keyword)
{
    std::vector<WellRefinementRecord> records;
    for (const auto& record : keyword) {
        WellRefinementRecord r;
        r.pattern = record.getItem("WELL").getTrimmedString(0);
        r.factor = { record.getItem("NX").get<int>(0),
                     record.getItem("NY").get<int>(0),
                     record.getItem("NZ").get<int>(0) };
        for (const char* ring : { "RING1", "RING2", "RING3" }) {
            const int w = record.getItem(ring).get<int>(0);
            if (w > 0) {
                r.rings.push_back(w);
            }
        }
        const auto layers = record.getItem("LAYERS").getTrimmedString(0);
        if (layers != "PERF" && layers != "ALL") {
            throw std::invalid_argument(
                fmt::format("WELLREF for '{}': LAYERS must be PERF or ALL, not '{}'.",
                            r.pattern, layers));
        }
        r.allLayers = (layers == "ALL");
        r.extraLayers = record.getItem("KLAYERS").get<int>(0);
        for (int f : r.factor) {
            if (f <= 0) {
                throw std::invalid_argument(
                    fmt::format("WELLREF for '{}': refinement factors must be positive.", r.pattern));
            }
        }
        if (r.rings.empty()) {
            throw std::invalid_argument(
                fmt::format("WELLREF for '{}': at least one ring width must be positive.", r.pattern));
        }
        records.push_back(std::move(r));
    }
    return records;
}

namespace {

struct Box
{
    std::array<int,3> lo{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() };
    std::array<int,3> hi{ std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min() };

    bool touches(const Box& o) const
    {
        for (int d = 0; d < 3; ++d) {
            if (hi[d] + 1 < o.lo[d] || o.hi[d] + 1 < lo[d]) {
                return false;
            }
        }
        return true;
    }
    void merge(const Box& o)
    {
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], o.lo[d]);
            hi[d] = std::max(hi[d], o.hi[d]);
        }
    }
    bool contains(const Box& o) const
    {
        for (int d = 0; d < 3; ++d) {
            if (o.lo[d] < lo[d] || o.hi[d] > hi[d]) {
                return false;
            }
        }
        return true;
    }
};

void mergeTouching(std::vector<Box>& boxes)
{
    bool merged = true;
    while (merged) {
        merged = false;
        for (std::size_t a = 0; a < boxes.size() && !merged; ++a) {
            for (std::size_t b = a + 1; b < boxes.size(); ++b) {
                if (boxes[a].touches(boxes[b])) {
                    boxes[a].merge(boxes[b]);
                    boxes.erase(boxes.begin() + b);
                    merged = true;
                    break;
                }
            }
        }
    }
}

} // anonymous namespace

std::vector<WellRefinementBox>
wellRefinementBoxes(const std::vector<WellRefinementRecord>& records,
                    const std::vector<std::pair<std::string, std::array<std::array<int,3>,2>>>& wellSeeds,
                    const std::array<int,3>& dims)
{
    std::vector<WellRefinementBox> out;
    int zoneIdx = 0;
    for (const auto& rec : records) {
        ++zoneIdx;
        std::vector<Box> seeds;
        for (const auto& [name, seed] : wellSeeds) {
            if (shmatch(rec.pattern, name)) {
                seeds.push_back({ seed[0], seed[1] });
            }
        }
        if (seeds.empty()) {
            OpmLog::warning(fmt::format("WELLREF: no completed well matches '{}'; ignored.", rec.pattern));
            continue;
        }

        const int nRings = static_cast<int>(rec.rings.size());
        std::vector<int> reach(nRings, 0);
        for (int r = nRings - 1; r >= 0; --r) {
            reach[r] = rec.rings[r] + ((r + 1 < nRings) ? reach[r + 1] : 0);
        }

        // Ring r (0 = outermost) is every seed dilated laterally by reach[r]
        // and, in K, by one layer per ring inside it, so that a nested ring is
        // interior to its parent in every direction; touching boxes merge.
        std::vector<std::vector<Box>> rings(nRings);
        for (int r = 0; r < nRings; ++r) {
            for (const auto& seed : seeds) {
                Box b;
                for (int d = 0; d < 2; ++d) {
                    b.lo[d] = std::max(0, seed.lo[d] - reach[r]);
                    b.hi[d] = std::min(dims[d] - 1, seed.hi[d] + reach[r]);
                }
                if (rec.allLayers) {
                    b.lo[2] = 0;
                    b.hi[2] = dims[2] - 1;
                }
                else {
                    const int kReach = rec.extraLayers + (nRings - 1 - r);
                    b.lo[2] = std::max(0, seed.lo[2] - kReach);
                    b.hi[2] = std::min(dims[2] - 1, seed.hi[2] + kReach);
                }
                rings[r].push_back(b);
            }
            mergeTouching(rings[r]);
        }

        struct Emitted { Box coarse; std::string name; std::array<int,3> unitsPerCoarse; };
        std::vector<Emitted> previous;
        for (int r = 0; r < nRings; ++r) {
            std::vector<Emitted> current;
            int boxIdx = 0;
            for (const auto& box : rings[r]) {
                ++boxIdx;
                WellRefinementBox wb;
                wb.name = fmt::format("WR{}R{}B{}", zoneIdx, r + 1, boxIdx);
                if (r == 0) {
                    for (int d = 0; d < 3; ++d) {
                        wb.lo[d] = box.lo[d] + 1;
                        wb.hi[d] = box.hi[d] + 1;
                        wb.nxyz[d] = (box.hi[d] - box.lo[d] + 1) * rec.factor[d];
                    }
                    current.push_back({ box, wb.name, rec.factor });
                    out.push_back(std::move(wb));
                    continue;
                }
                const Emitted* parent = nullptr;
                for (const auto& p : previous) {
                    if (p.coarse.contains(box)) {
                        parent = &p;
                        break;
                    }
                }
                if (parent == nullptr) {
                    throw std::logic_error("WELLREF: an inner ring is not contained in any outer ring box");
                }
                // The grid builder needs a nested block strictly interior to its
                // parent; clamping at the grid edge can take that away.
                Box inner = box;
                bool shrunk = false;
                for (int d = 0; d < 3; ++d) {
                    if (inner.lo[d] <= parent->coarse.lo[d]) { inner.lo[d] = parent->coarse.lo[d] + 1; shrunk = true; }
                    if (inner.hi[d] >= parent->coarse.hi[d]) { inner.hi[d] = parent->coarse.hi[d] - 1; shrunk = true; }
                }
                if (inner.lo[0] > inner.hi[0] || inner.lo[1] > inner.hi[1] || inner.lo[2] > inner.hi[2]) {
                    OpmLog::warning(fmt::format("WELLREF '{}': ring {} would touch its parent's "
                                                "boundary at the grid edge and is dropped.",
                                                rec.pattern, r + 1));
                    continue;
                }
                if (shrunk) {
                    OpmLog::warning(fmt::format("WELLREF '{}': ring {} shrunk to stay inside its "
                                                "parent at the grid edge.", rec.pattern, r + 1));
                }
                wb.parent = parent->name;
                std::array<int,3> units{};
                for (int d = 0; d < 3; ++d) {
                    // Parent-local, 1-based, in the parent's refined cells.
                    wb.lo[d] = (inner.lo[d] - parent->coarse.lo[d]) * parent->unitsPerCoarse[d] + 1;
                    wb.hi[d] = (inner.hi[d] + 1 - parent->coarse.lo[d]) * parent->unitsPerCoarse[d];
                    wb.nxyz[d] = (wb.hi[d] - wb.lo[d] + 1) * rec.factor[d];
                    units[d] = parent->unitsPerCoarse[d] * rec.factor[d];
                }
                current.push_back({ inner, wb.name, units });
                out.push_back(std::move(wb));
            }
            previous = std::move(current);
        }
    }
    return out;
}

std::vector<WellRefinementBox> wellRefinementBoxes(const Deck& deck)
{
    std::vector<WellRefinementRecord> records;
    for (const auto* kw : deck.getKeywordList("WELLREF")) {
        auto r = readWellRefinement(*kw);
        records.insert(records.end(), r.begin(), r.end());
    }
    if (records.empty()) {
        return {};
    }
    if (!deck.hasKeyword("DIMENS")) {
        throw std::invalid_argument("WELLREF needs DIMENS.");
    }
    const auto& dimens = deck["DIMENS"].back().getRecord(0);
    const std::array<int,3> dims{ dimens.getItem("NX").get<int>(0),
                                  dimens.getItem("NY").get<int>(0),
                                  dimens.getItem("NZ").get<int>(0) };

    // Well heads, then every COMPDAT record's cells; defaulted I/J take the head.
    std::map<std::string, std::array<int,2>> heads;
    for (const auto* kw : deck.getKeywordList("WELSPECS")) {
        for (const auto& rec : *kw) {
            heads[rec.getItem("WELL").getTrimmedString(0)] =
                { rec.getItem("HEAD_I").get<int>(0) - 1, rec.getItem("HEAD_J").get<int>(0) - 1 };
        }
    }
    std::map<std::string, std::array<std::array<int,3>,2>> seeds;
    for (const auto* kw : deck.getKeywordList("COMPDAT")) {
        for (const auto& rec : *kw) {
            const auto wellPattern = rec.getItem("WELL").getTrimmedString(0);
            for (const auto& [name, head] : heads) {
                if (!shmatch(wellPattern, name)) {
                    continue;
                }
                const auto& iItem = rec.getItem("I");
                const auto& jItem = rec.getItem("J");
                const int i = (iItem.defaultApplied(0) || iItem.get<int>(0) <= 0) ? head[0] : iItem.get<int>(0) - 1;
                const int j = (jItem.defaultApplied(0) || jItem.get<int>(0) <= 0) ? head[1] : jItem.get<int>(0) - 1;
                const int k1 = rec.getItem("K1").get<int>(0) - 1;
                const int k2 = rec.getItem("K2").get<int>(0) - 1;
                auto it = seeds.find(name);
                if (it == seeds.end()) {
                    it = seeds.emplace(name, std::array<std::array<int,3>,2>{
                        std::array<int,3>{ i, j, k1 }, std::array<int,3>{ i, j, k2 } }).first;
                }
                else {
                    auto& s = it->second;
                    s[0] = { std::min(s[0][0], i), std::min(s[0][1], j), std::min(s[0][2], k1) };
                    s[1] = { std::max(s[1][0], i), std::max(s[1][1], j), std::max(s[1][2], k2) };
                }
            }
        }
    }
    for (const auto* kwName : { "COMPTRAJ", "WELTRAJ" }) {
        if (deck.hasKeyword(kwName)) {
            OpmLog::warning(fmt::format("WELLREF: wells completed by {} are not used to place "
                                        "refinement; only COMPDAT wells are.", kwName));
        }
    }
    std::vector<std::pair<std::string, std::array<std::array<int,3>,2>>> seedList(seeds.begin(), seeds.end());
    return wellRefinementBoxes(records, seedList, dims);
}

} // namespace Opm
