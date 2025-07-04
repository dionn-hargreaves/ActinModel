// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.4 (2016/08/29)

#pragma once

#include <Mathematics/GteMesh.h>
#include <Mathematics/GteRectangle.h>
#include <Mathematics/GteVector3.h>
#include <memory>

namespace gte
{

template <typename Real>
class RectangleMesh : public Mesh<Real>
{
public:
    // Create a mesh that tessellates a rectangle.
    RectangleMesh(MeshDescription const& description, Rectangle<3, Real> const& rectangle);

    // Member access.
    inline Rectangle<3, Real> const& GetRectangle() const;

protected:
    void InitializeTCoords();
    void InitializePositions();
    void InitializeNormals();
    void InitializeFrame();

    Rectangle<3, Real> mRectangle;

    // If the client does not request texture coordinates, they will be
    // computed internally for use in evaluation of the surface geometry.
    std::vector<Vector2<Real>> mDefaultTCoords;
};


template <typename Real>
RectangleMesh<Real>::RectangleMesh(MeshDescription const& description,
    Rectangle<3, Real> const& rectangle)
    :
    Mesh<Real>(description, { MeshTopology::RECTANGLE }),
    mRectangle(rectangle)
{
    if (!this->mDescription.constructed)
    {
        // The logger system will report these errors in the Mesh constructor.
        return;
    }

    if (!this->mTCoords)
    {
        mDefaultTCoords.resize(this->mDescription.numVertices);
        this->mTCoords = mDefaultTCoords.data();
        this->mTCoordStride = sizeof(Vector2<Real>);

        this->mDescription.allowUpdateFrame = this->mDescription.wantDynamicTangentSpaceUpdate;
        if (this->mDescription.allowUpdateFrame)
        {
            if (!this->mDescription.hasTangentSpaceVectors)
            {
                this->mDescription.allowUpdateFrame = false;
            }

            if (!this->mNormals)
            {
                this->mDescription.allowUpdateFrame = false;
            }
        }
    }

    this->ComputeIndices();
    InitializeTCoords();
    InitializePositions();
    if (this->mDescription.allowUpdateFrame)
    {
        InitializeFrame();
    }
    else if (this->mNormals)
    {
        InitializeNormals();
    }
}

template <typename Real>
inline Rectangle<3, Real> const& RectangleMesh<Real>::GetRectangle() const
{
    return mRectangle;
}

template <typename Real>
void RectangleMesh<Real>::InitializeTCoords()
{
    Vector2<Real> tcoord;
    for (uint32_t r = 0, i = 0; r < this->mDescription.numRows; ++r)
    {
        tcoord[1] = (Real)r / (Real)(this->mDescription.numRows - 1);
        for (uint32_t c = 0; c < this->mDescription.numCols; ++c, ++i)
        {
            tcoord[0] = (Real)c / (Real)(this->mDescription.numCols - 1);
            this->TCoord(i) = tcoord;
        }
    }
}

template <typename Real>
void RectangleMesh<Real>::InitializePositions()
{
    for (uint32_t r = 0, i = 0; r < this->mDescription.numRows; ++r)
    {
        for (uint32_t c = 0; c < this->mDescription.numCols; ++c, ++i)
        {
            Vector2<Real> tcoord = this->TCoord(i);
            Real w0 = ((Real)2 * tcoord[0] - (Real)1) * mRectangle.extent[0];
            Real w1 = ((Real)2 * tcoord[1] - (Real)1) * mRectangle.extent[1];
            this->Position(i) = mRectangle.center + w0 * mRectangle.axis[0] + w1 * mRectangle.axis[1];
        }
    }
}

template <typename Real>
void RectangleMesh<Real>::InitializeNormals()
{
    Vector3<Real> normal = UnitCross(mRectangle.axis[0], mRectangle.axis[1]);
    for (uint32_t i = 0; i < this->mDescription.numVertices; ++i)
    {
        this->Normal(i) = normal;
    }
}

template <typename Real>
void RectangleMesh<Real>::InitializeFrame()
{
    Vector3<Real> normal = UnitCross(mRectangle.axis[0], mRectangle.axis[1]);
    Vector3<Real> tangent{ (Real)1, (Real)0, (Real)0 };
    Vector3<Real> bitangent{ (Real)0, (Real)1, (Real)0 };  // = Cross(normal,tangent)
    for (uint32_t i = 0; i < this->mDescription.numVertices; ++i)
    {
        if (this->mNormals)
        {
            this->Normal(i) = normal;
        }

        if (this->mTangents)
        {
            this->Tangent(i) = tangent;
        }

        if (this->mBitangents)
        {
            this->Bitangent(i) = bitangent;
        }

        if (this->mDPDUs)
        {
            this->DPDU(i) = tangent;
        }

        if (this->mDPDVs)
        {
            this->DPDV(i) = bitangent;
        }
    }
}

}
