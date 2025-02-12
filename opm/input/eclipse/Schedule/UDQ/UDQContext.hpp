/*
  Copyright 2019 Equinor ASA.

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

#ifndef UDQ_CONTEXT_HPP
#define UDQ_CONTEXT_HPP

#include <opm/input/eclipse/EclipseState/Grid/RegionSetMatcher.hpp>
#include <opm/input/eclipse/Schedule/MSW/SegmentMatcher.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {

    class GroupOrder;
    class SegmentSet;
    class SummaryState;
    class UDQFunctionTable;
    class UDQSet;
    class UDQState;
    class UDT;
    class WellMatcher;

} // namespace Opm

namespace Opm {

    class UDQContext
    {
    public:
        struct MatcherFactories
        {
            std::function<std::unique_ptr<SegmentMatcher>()> segments{};
            std::function<std::unique_ptr<RegionSetMatcher>()> regions{};
        };

        UDQContext(const UDQFunctionTable& udqft,
                   const WellMatcher&      wm,
                   const GroupOrder&       go,
                   const std::unordered_map<std::string, UDT>& tables,
                   MatcherFactories        create_matchers,
                   SummaryState&           summary_state,
                   UDQState&               udq_state);

        std::optional<double> get(const std::string& key) const;

        std::optional<double>
        get_well_var(const std::string& well, const std::string& var) const;

        std::optional<double>
        get_group_var(const std::string& group, const std::string& var) const;

        std::optional<double>
        get_segment_var(const std::string& well,
                        const std::string& var,
                        std::size_t segment) const;

        std::optional<double>
        get_region_var(const std::string& regSet,
                       const std::string& var,
                       std::size_t region) const;

        const UDT& get_udt(const std::string& name) const;

        void add(const std::string& key, double value);
        void update_assign(const std::string& keyword, const UDQSet& udq_result);
        void update_define(std::size_t report_step,
                           const std::string& keyword,
                           const UDQSet& udq_result);

        const UDQFunctionTable& function_table() const;

        const std::vector<std::string>& wells() const;
        std::vector<std::string> wells(const std::string& pattern) const;
        std::vector<std::string> nonFieldGroups() const;
        std::vector<std::string> groups(const std::string& pattern) const;
        SegmentSet segments() const;
        SegmentSet segments(const std::vector<std::string>& set_descriptor) const;

        RegionSetMatchResult regions() const;
        RegionSetMatchResult regions(const std::string&              vector_name,
                                     const std::vector<std::string>& set_descriptor) const;

    private:
        struct Matchers
        {
            std::unique_ptr<SegmentMatcher> segments{};
            std::unique_ptr<RegionSetMatcher> regions{};
        };

        const UDQFunctionTable& udqft;
        const WellMatcher& well_matcher;
        const GroupOrder& group_order_;
        const std::unordered_map<std::string, UDT>& udt;

        SummaryState& summary_state;
        UDQState& udq_state;

        MatcherFactories create_matchers_{};
        mutable Matchers matchers_{};

        //std::unordered_map<std::string, UDQSet> udq_results;
        std::unordered_map<std::string, double> values;

        void ensure_segment_matcher_exists() const;
        void ensure_region_matcher_exists() const;
    };

} // namespace Opm

#endif // UDQ_CONTEXT_HPP
