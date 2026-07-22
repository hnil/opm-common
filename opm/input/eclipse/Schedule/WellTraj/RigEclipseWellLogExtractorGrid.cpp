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

#include "RigEclipseWellLogExtractorGrid.hpp"

#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>
#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>

#include <external/resinsight/CommonCode/cvfStructGrid.h>

#include <external/resinsight/LibGeometry/cvfBoundingBox.h>

#include <array>
#include <cstddef>
#include <map>
#include <utility>

namespace external {

RigEclipseWellLogExtractorGrid::RigEclipseWellLogExtractorGrid( const RigWellPath* wellpath,
                                                               std::vector<std::array<cvf::Vec3d, 8>> cellCorners,
                                                               cvf::ref<cvf::BoundingBoxTree>& cellSearchTree )
    : RigWellLogExtractor( wellpath, "" )
      ,m_cellCorners( std::move( cellCorners ) )
      ,m_cellSearchTree( cellSearchTree )
{
    calculateIntersection();
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp
void RigEclipseWellLogExtractorGrid::calculateIntersection()
{
    std::map<RigMDCellIdxEnterLeaveKey, HexIntersectionInfo> uniqueIntersections;

    if ( m_wellPathGeometry->wellPathPoints().empty() ) return;

    buildCellSearchTree();

    for ( std::size_t wpp = 0; wpp < m_wellPathGeometry->wellPathPoints().size() - 1; ++wpp )
    {
        std::vector<HexIntersectionInfo> intersections;
        cvf::Vec3d                       p1 = m_wellPathGeometry->wellPathPoints()[wpp];
        cvf::Vec3d                       p2 = m_wellPathGeometry->wellPathPoints()[wpp + 1];

        cvf::BoundingBox bb;

        bb.add( p1 );
        bb.add( p2 );

        std::vector<std::size_t> closeCellIndices = findCloseCellIndices( bb );

        for ( const auto& globalCellIndex : closeCellIndices )
        {
            cvf::Vec3d    hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
            this->hexCornersOpmToResinsight( hexCorners, globalCellIndex);
            RigHexIntersectionTools::lineHexCellIntersection( p1, p2, hexCorners, globalCellIndex, &intersections );
        }

        // Now, with all the intersections of this piece of line, we need to
        // sort them in order, and set the measured depth and corresponding cell index

        // Inserting the intersections in this map will remove identical intersections
        // and sort them according to MD, CellIdx, Leave/enter

        double md1 = m_wellPathGeometry->measuredDepths()[wpp];
        double md2 = m_wellPathGeometry->measuredDepths()[wpp + 1];

        insertIntersectionsInMap( intersections, p1, md1, p2, md2, &uniqueIntersections );
    }

    this->populateReturnArrays( uniqueIntersections );
}

cvf::Vec3d RigEclipseWellLogExtractorGrid::calculateLengthInCell( std::size_t       cellIndex,
                                                                 const cvf::Vec3d& startPoint,
                                                                 const cvf::Vec3d& endPoint ) const
{
    cvf::Vec3d hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
    this->hexCornersOpmToResinsight( hexCorners, cellIndex );

    std::array<cvf::Vec3d, 8> hexCorners2;
    for (std::size_t i = 0; i < 8; i++)
           hexCorners2[i] =  hexCorners[i];

    return this->calculateLengthInCell( hexCorners2, startPoint, endPoint );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigWellPathIntersectionTools.cpp
cvf::Vec3d RigEclipseWellLogExtractorGrid::calculateLengthInCell( const std::array<cvf::Vec3d, 8>& hexCorners,
                                                                 const cvf::Vec3d&                startPoint,
                                                                 const cvf::Vec3d&                endPoint ) const
{

    cvf::Vec3d vec = endPoint - startPoint;
    cvf::Vec3d iAxisDirection;
    cvf::Vec3d jAxisDirection;
    cvf::Vec3d kAxisDirection;

    this->findCellLocalXYZ( hexCorners, iAxisDirection, jAxisDirection, kAxisDirection );

    cvf::Mat3d localCellCoordinateSystem( iAxisDirection.x(),
                                          jAxisDirection.x(),
                                          kAxisDirection.x(),
                                          iAxisDirection.y(),
                                          jAxisDirection.y(),
                                          kAxisDirection.y(),
                                          iAxisDirection.z(),
                                          jAxisDirection.z(),
                                          kAxisDirection.z() );

    auto signedVector = vec.getTransformedVector( localCellCoordinateSystem.getInverted() );

    return { std::fabs( signedVector.x() ), std::fabs( signedVector.y() ), std::fabs( signedVector.z() ) };
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigCellGeometryTools.cpp
void RigEclipseWellLogExtractorGrid::findCellLocalXYZ( const std::array<cvf::Vec3d, 8>& hexCorners,
                                                      cvf::Vec3d&                      localXdirection,
                                                      cvf::Vec3d&                      localYdirection,
                                                      cvf::Vec3d&                      localZdirection ) const
{

    cvf::Vec3d faceCenterNegI = cvf::GeometryTools::computeFaceCenter( hexCorners[0],
                                                                       hexCorners[4],
                                                                       hexCorners[7],
                                                                       hexCorners[3] );

    cvf::Vec3d faceCenterPosI = cvf::GeometryTools::computeFaceCenter( hexCorners[1],
                                                                       hexCorners[2],
                                                                       hexCorners[6],
                                                                       hexCorners[5] );

    cvf::Vec3d faceCenterNegJ = cvf::GeometryTools::computeFaceCenter( hexCorners[0],
                                                                       hexCorners[1],
                                                                       hexCorners[5],
                                                                       hexCorners[4] );

    cvf::Vec3d faceCenterPosJ = cvf::GeometryTools::computeFaceCenter( hexCorners[3],
                                                                       hexCorners[7],
                                                                       hexCorners[6],
                                                                       hexCorners[2] );

    cvf::Vec3d faceCenterCenterVectorI = faceCenterPosI - faceCenterNegI;
    cvf::Vec3d faceCenterCenterVectorJ = faceCenterPosJ - faceCenterNegJ;

    localZdirection.cross( faceCenterCenterVectorI, faceCenterCenterVectorJ );
    localZdirection.normalize();

    cvf::Vec3d crossPoductJZ;
    crossPoductJZ.cross( faceCenterCenterVectorJ, localZdirection );
    localXdirection = faceCenterCenterVectorI + crossPoductJZ;
    localXdirection.normalize();

    cvf::Vec3d crossPoductIZ;
    crossPoductIZ.cross( faceCenterCenterVectorI, localZdirection );
    localYdirection = faceCenterCenterVectorJ - crossPoductIZ;
    localYdirection.normalize();
}

// The eight corner points of each cell are supplied in OPM/ECL getCornerPos()
// order (l = 0..7); convert to the ResInsight hex numbering used by the
// intersection routines.
void RigEclipseWellLogExtractorGrid::hexCornersOpmToResinsight( cvf::Vec3d hexCorners[8], std::size_t cellIndex ) const
{
    static const std::array<std::size_t, 8> opm2resinsight = {0, 1, 3, 2, 4, 5, 7, 6};

    const auto& corners = m_cellCorners[cellIndex];
    for (std::size_t l = 0; l < 8; l++) {
        hexCorners[opm2resinsight[l]] = corners[l];
    }
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractorGrid::buildCellSearchTree()
{
    if (m_cellSearchTree.isNull()) {

        const std::size_t cellCount = m_cellCorners.size();

        std::vector<std::size_t>      cellIndicesForBoundingBoxes;
        std::vector<cvf::BoundingBox> cellBoundingBoxes;

        cellIndicesForBoundingBoxes.reserve( cellCount );
        cellBoundingBoxes.reserve( cellCount );

        for (std::size_t cIdx = 0; cIdx < cellCount; ++cIdx) {
            const auto& corners = m_cellCorners[cIdx];
            cvf::BoundingBox cellBB;
            for (std::size_t l = 0; l < 8; l++) {
                cellBB.add( corners[l] );
            }

            if (cellBB.isValid()) {
                cellIndicesForBoundingBoxes.emplace_back(cIdx);
                cellBoundingBoxes.emplace_back(cellBB);
            }
        }

        cellIndicesForBoundingBoxes.shrink_to_fit();
        cellBoundingBoxes.shrink_to_fit();

        m_cellSearchTree = new cvf::BoundingBoxTree;
        m_cellSearchTree->buildTreeFromBoundingBoxes( cellBoundingBoxes, &cellIndicesForBoundingBoxes );
    }
}

// From ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractorGrid::findIntersectingCells( const cvf::BoundingBox& inputBB, std::vector<std::size_t>* cellIndices ) const
{
    CVF_ASSERT( m_cellSearchTree.notNull() );

    m_cellSearchTree->findIntersections( inputBB, cellIndices );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp
std::vector<std::size_t> RigEclipseWellLogExtractorGrid::findCloseCellIndices( const cvf::BoundingBox& bb )
{
    std::vector<std::size_t> closeCells;
    this->findIntersectingCells( bb, &closeCells );
    return closeCells;
}

cvf::ref<cvf::BoundingBoxTree> RigEclipseWellLogExtractorGrid::getCellSearchTree()
{
    return m_cellSearchTree;
}

} //namespace external
