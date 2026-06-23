/*
  Copyright 2026 Equinor.

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

#pragma once

#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/LibGeometry/cvfBoundingBoxTree.h>

#include <array>
#include <cstddef>
#include <vector>

namespace external {

class RigWellPath;

namespace cvf
{
class BoundingBox;
}

//==================================================================================================
/// Trajectory/grid intersection driven by *explicit cell geometry* rather than
/// an Opm::EclipseGrid.
///
/// This is the geometry-only sibling of RigEclipseWellLogExtractor: it builds
/// the same cvf AABB search tree and runs the same line/hex intersection, but
/// the cells are supplied as a flat list of corner points. That lets a WELTRAJ
/// trajectory be intersected against an arbitrary grid — in particular a
/// locally refined (LGR) leaf grid — as a post-process, instead of only the
/// coarse EclipseGrid used by RigEclipseWellLogExtractor / loadCOMPTRAJ.
///
/// `cellCorners[c]` holds the eight corner points of cell `c` in the OPM/ECL
/// getCornerPos() order (local corner index l = 0..7). The OPM->ResInsight hex
/// permutation {0,1,3,2,4,5,7,6} is applied internally, exactly as in
/// RigEclipseWellLogExtractor, so callers pass corners in the natural OPM order.
/// The intersection results reference cells by their index into `cellCorners`.
//==================================================================================================
class RigEclipseWellLogExtractorGrid : public RigWellLogExtractor
{
public:
    RigEclipseWellLogExtractorGrid( const RigWellPath* wellpath,
                                    std::vector<std::array<cvf::Vec3d, 8>> cellCorners,
                                    cvf::ref<cvf::BoundingBoxTree>& cellSearchTree );

    cvf::ref<cvf::BoundingBoxTree> getCellSearchTree();

private:
    void                     calculateIntersection();
    std::vector<std::size_t> findCloseCellIndices( const cvf::BoundingBox& bb );

    cvf::Vec3d calculateLengthInCell( std::size_t       cellIndex,
                                      const cvf::Vec3d& startPoint,
                                      const cvf::Vec3d& endPoint ) const override;

    cvf::Vec3d calculateLengthInCell( const std::array<cvf::Vec3d, 8>& hexCorners,
                                      const cvf::Vec3d&                startPoint,
                                      const cvf::Vec3d&                endPoint ) const;

    void hexCornersOpmToResinsight( cvf::Vec3d hexCorners[8], std::size_t cellIndex ) const;

    void findCellLocalXYZ( const std::array<cvf::Vec3d, 8>& hexCorners,
                           cvf::Vec3d&                      localXdirection,
                           cvf::Vec3d&                      localYdirection,
                           cvf::Vec3d&                      localZdirection ) const;

    void buildCellSearchTree();
    void findIntersectingCells( const cvf::BoundingBox& inputBB, std::vector<std::size_t>* cellIndices ) const;

    std::vector<std::array<cvf::Vec3d, 8>> m_cellCorners;
    cvf::ref<cvf::BoundingBoxTree>         m_cellSearchTree;
};

} // namespace external
