/*
  Copyright 2013 Statoil ASA.

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

#ifndef CONNECTIONSET_HPP_
#define CONNECTIONSET_HPP_

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Opm {
    class ActiveGridCells;
    class DeckRecord;
    class EclipseGrid;
    class ErrorGuard;
    class FieldPropsManager;
    class KeywordLocation;
    class ParseContext;
    class ScheduleGrid;
    class WDFAC;
    struct WellTrajInfo;
} // namespace Opm

namespace Opm {

    class WellConnections
    {
    public:
        using const_iterator = std::vector<Connection>::const_iterator;

        WellConnections() = default;
        WellConnections(const Connection::Order ordering, const int headI, const int headJ);
        WellConnections(const Connection::Order ordering, const int headI, const int headJ,
                        const std::vector<Connection>& connections);

        static WellConnections serializationTestObject();

        // cppcheck-suppress noExplicitConstructor
        template <class Grid>
        WellConnections(const WellConnections& src, const Grid& grid)
            : m_ordering(src.ordering())
            , headI     (src.headI)
            , headJ     (src.headJ)
        {
            for (const auto& c : src) {
                if (grid.isCellActive(c.getI(), c.getJ(), c.getK())) {
                    this->add(c);
                }
            }
        }

        void add(const Connection& conn)
        {
            this->m_connections.push_back(conn);
        }

        void addConnection(const int i, const int j, const int k,
                           const std::size_t global_index,
                           const Connection::State state,
                           const double depth,
                           const Connection::CTFProperties& ctf_props,
                           const int satTableId,
                           const Connection::Direction direction = Connection::Direction::Z,
                           const Connection::CTFKind ctf_kind = Connection::CTFKind::DeckValue,
                           const std::size_t seqIndex = 0,
                           int lgr_grid_number = 0,
                           const bool defaultSatTabId = true);

        // Out-parameter 'requested_open_complnums' is cleared and then set to
        // the complnum of every existing connection this record opens (assigns
        // OPEN), for raising REQUEST_OPEN_COMPLETION events.  Likewise,
        // 'requested_shut_complnums' is cleared and then set to the complnum of
        // every existing connection this record shuts (assigns a non-OPEN
        // state), for clearing any pending REQUEST_OPEN_COMPLETION events.
        void loadCOMPDAT(const DeckRecord&      record,
                         const std::string&     wname,
                         const WDFAC&           wdfac,
                         const ScheduleGrid&    grid,
                         const KeywordLocation& location,
                         const ParseContext&    parseContext,
                         ErrorGuard&            errors,
                         std::vector<int>&      requested_open_complnums,
                         std::vector<int>&      requested_shut_complnums);

        void loadCOMPDATL(const DeckRecord&      record,
                          const std::string&     wname,
                          const WDFAC&           wdfac,
                          const ScheduleGrid&    grid,
                          const KeywordLocation& location,
                          const ParseContext&    parseContext,
                          ErrorGuard&            errors,
                          std::vector<int>&      requested_open_complnums,
                          std::vector<int>&      requested_shut_complnums);

        void loadCOMPTRAJ(const DeckRecord&      record,
                          const std::string&     wname,
                          const ScheduleGrid&    grid,
                          const KeywordLocation& location,
                          WellTrajInfo&          wellTraj);

        void loadWELTRAJ(const DeckRecord&      record,
                         const std::string&     wname,
                         const ScheduleGrid&    grid,
                         const KeywordLocation& location);

        void applyDFactorCorrelation(const ScheduleGrid& grid,
                                     const WDFAC&        wdfac);

        int getHeadI() const;
        int getHeadJ() const;
        const std::vector<double>& getMD() const;

        /// Whether this connection set was generated from a well trajectory
        /// (WELTRAJ/COMPTRAJ) rather than from explicit COMPDAT cells.
        bool hasTrajectory() const;

        /// Per-cell geometry/properties needed to (re)build a trajectory
        /// connection against an arbitrary grid. Supplied by the caller of
        /// recomputeTrajectoryConnections() for each intersected cell.
        struct TrajectoryCell {
            std::array<int, 3>    ijk{};            //!< logical cartesian index recorded on the connection
            std::size_t           global_index{};   //!< grid cell index recorded on the connection
            double                depth{};          //!< cell centre depth
            std::array<double, 3> perm{};           //!< permx, permy, permz
            std::array<double, 3> dimensions{};     //!< cell extent dx, dy, dz
            double                ntg{1.0};
            int                   satnum{0};
        };

        /// Rebuild the trajectory connections of this well by intersecting the
        /// retained WELTRAJ trajectory against an explicitly supplied grid
        /// geometry -- e.g. a locally refined (LGR) leaf grid -- instead of the
        /// coarse EclipseGrid used at parse time (loadCOMPTRAJ).
        ///
        /// \param[in] cellCorners  per-cell eight corner points (OPM/ECL
        ///            getCornerPos order), indexed by the cell index that
        ///            cellInfo() is queried with.
        /// \param[in] cellInfo  returns the properties of the intersected cell
        ///            with the given index, or std::nullopt if it must be
        ///            skipped (e.g. inactive). The recorded connection uses the
        ///            ijk/global_index from this callback.
        ///
        /// Does nothing if this well has no retained trajectory. Reuses the
        /// same CTF/Kh computation as loadCOMPTRAJ.
        void recomputeTrajectoryConnections
            (const std::vector<std::array<std::array<double,3>, 8>>&                 cellCorners,
             const std::function<std::optional<TrajectoryCell>(std::size_t)>&        cellInfo);

        std::size_t size() const;
        bool empty() const;
        std::size_t num_open() const;
        const Connection& operator[](std::size_t index) const;
        const Connection& get(std::size_t index) const;
        const Connection& getFromIJK(const int i, const int j, const int k) const;
        const Connection& getFromGlobalIndex(std::size_t global_index) const;
        const Connection& lowest() const;
        Connection& getFromIJK(const int i, const int j, const int k);
        Connection* maybeGetFromGlobalIndex(const std::size_t global_index);
        bool hasGlobalIndex(std::size_t global_index) const;
        double segment_perf_length(int segment) const;

        const_iterator begin() const { return this->m_connections.begin(); }
        const_iterator end() const { return this->m_connections.end(); }
        auto begin() { return this->m_connections.begin(); }
        auto end() { return this->m_connections.end(); }
        bool allConnectionsShut() const;
        /// Order connections irrespective of input order.
        /// The algorithm used is the following:
        ///     1. The connection nearest to the given (well_i, well_j)
        ///        coordinates in terms of the connection's (i, j) is chosen
        ///        to be the first connection. If non-unique, choose one with
        ///        lowest z-depth (shallowest).
        ///     2. Choose next connection to be nearest to current in (i, j) sense.
        ///        If non-unique choose closest in z-depth (not logical cartesian k).
        ///
        /// \param[in] well_i  logical cartesian i-coordinate of well head
        /// \param[in] well_j  logical cartesian j-coordinate of well head
        /// \param[in] grid    EclipseGrid object, used for cell depths
        void order();

        bool operator==( const WellConnections& ) const;
        bool operator!=( const WellConnections& ) const;

        Connection::Order ordering() const { return this->m_ordering; }
        std::vector<const Connection *> output(const EclipseGrid& grid) const;

        /// Activate or reactivate WELPI scaling for this connection set.
        ///
        /// Following this call, any WELPI-based scaling will apply to all
        /// connections whose properties are not reset in COMPDAT.
        ///
        /// Returns whether or not this call to prepareWellPIScaling() is
        /// a state change (e.g., no WELPI to active WELPI or WELPI for
        /// some connections to WELPI for all connections).
        bool prepareWellPIScaling();

        /// Scale pertinent connections' CF value by supplied value.  Scaling
        /// factor typically derived from 'WELPI' input keyword and a dynamic
        /// productivity index calculation.  Applicability array specifies
        /// whether or not a particular connection is exempt from scaling.
        /// Empty array means "apply scaling to all eligible connections".
        /// This array is updated on return (entries set to 'false' if
        /// corresponding connection is not eligible).
        void applyWellPIScaling(const double       scaleFactor,
                                std::vector<bool>& scalingApplicable);

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->m_ordering);
            serializer(this->headI);
            serializer(this->headJ);
            serializer(this->m_connections);
            serializer(this->coord);
            serializer(this->md);
            serializer(this->m_traj_perfs);
        }

    private:
        Connection::Order m_ordering { Connection::Order::TRACK };
        int headI{0};
        int headJ{0};
        std::vector<Connection> m_connections{};

        std::array<std::vector<double>, 3> coord{};
        std::vector<double> md{};

        /// Per-COMPTRAJ-record parameters retained at parse time so the
        /// trajectory connections can be recomputed against a different (e.g.
        /// refined) grid. Mirrors the items consumed by loadCOMPTRAJ.
        struct TrajPerf {
            double perf_top{};
            double perf_bot{};
            double rw{};
            double skin_factor{};
            double d_factor{};
            double user_Kh{-1.0};        //!< explicit Kh, or < 0 if defaulted
            double user_CF{-1.0};        //!< explicit CF, or < 0 if defaulted
            int    sat_table_id{-1};     //!< explicit SATNUM table, or < 0 if defaulted
            bool   default_sat_table{true};
            Connection::State state{Connection::State::OPEN};

            template <class Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(perf_top); serializer(perf_bot);
                serializer(rw); serializer(skin_factor); serializer(d_factor);
                serializer(user_Kh); serializer(user_CF);
                serializer(sat_table_id); serializer(default_sat_table);
                serializer(state);
            }
            bool operator==(const TrajPerf&) const = default;
        };
        std::vector<TrajPerf> m_traj_perfs{};

        /// Shared per-cell CTF computation + connection add/update, used by both
        /// loadCOMPTRAJ (EclipseGrid) and recomputeTrajectoryConnections (grid).
        void addOrUpdateTrajectoryConnection(const TrajPerf&            rec,
                                             const TrajectoryCell&      cell,
                                             const std::array<double,3>& connection_vector);

        void addConnection(const int i, const int j, const int k,
                           const std::size_t global_index,
                           const int complnum,
                           const Connection::State state,
                           const double depth,
                           const Connection::CTFProperties& ctf_props,
                           const int satTableId,
                           const Connection::Direction direction,
                           const Connection::CTFKind ctf_kind,
                           const std::size_t seqIndex,
                           int lgr_grid_number,
                           const bool defaultSatTabId);

        std::size_t findClosestConnection(int oi, int oj, double oz, std::size_t start_pos);
        void orderTRACK();
        void orderMSW();
        void orderDEPTH();

        void loadCOMPDATX(const DeckRecord&                 record,
                          const std::string&                wname,
                          const WDFAC&                      wdfac,
                          const ScheduleGrid&               grid,
                          const KeywordLocation&            location,
                          const std::optional<std::string>& lgr_label,
                          const ParseContext&               parseContext,
                          ErrorGuard&                       errors,
                          std::vector<int>&                 requested_open_complnums,
                          std::vector<int>&                 requested_shut_complnums);
    };

    std::optional<int>
    getCompletionNumberFromGlobalConnectionIndex(const WellConnections& connections,
                                                 const std::size_t      global_index);
} // namespace Opm

#endif // CONNECTIONSET_HPP_
