// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteExtremalQuery3.h>
#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteVETManifoldMesh.h>
#include <LowLevel/GteLogger.h>
#include <LowLevel/GteRangeIteration.h>
#include <algorithm>
#include <stack>
#include <queue>

namespace gte
{

template <typename Real>
class ExtremalQuery3BSP : public ExtremalQuery3<Real>
{
public:
    // Construction.
    ExtremalQuery3BSP(Polyhedron3<Real> const& polytope);

    // Disallow copying and assignment.
    ExtremalQuery3BSP(ExtremalQuery3BSP const&) = delete;
    ExtremalQuery3BSP& operator=(ExtremalQuery3BSP const&) = delete;

    // Compute the extreme vertices in the specified direction and return the
    // indices of the vertices in the polyhedron vertex array.
    virtual void GetExtremeVertices(Vector3<Real> const& direction,
        int& positiveDirection, int& negativeDirection) override;

    // Tree statistics.
    inline int GetNumNodes() const;
    inline int GetTreeDepth() const;

private:
    class SphericalArc
    {
    public:
        // Construction.
        SphericalArc();

        // The arcs are stored in a multiset ordered by increasing separation.
        // The multiset will be traversed in reverse order.  This heuristic is
        // designed to create BSP trees whose top-most nodes can eliminate as
        // many arcs as possible during an extremal query.
        bool operator<(SphericalArc const& arc) const;

        // Indices N[] into the face normal array for the endpoints of the arc.
        int nIndex[2];

        // The number of arcs in the path from normal N[0] to normal N[1].
        // For spherical polygon edges, the number is 1.  The number is 2 or
        // larger for bisector arcs of the spherical polygon.
        int separation;

        // The normal is Cross(FaceNormal[N[0]],FaceNormal[N[1]]).
        Vector3<Real> normal;

        // Indices into the vertex array for the extremal points for the
        // two regions sharing the arc.  As the arc is traversed from normal
        // N[0] to normal N[1], PosVertex is the index for the extreme vertex
        // to the left of the arc and NegVertex is the index for the extreme
        // vertex to the right of the arc.
        int posVertex, negVertex;

        // Support for BSP trees stored as contiguous nodes in an array.
        int posChild, negChild;
    };

    typedef VETManifoldMesh::Triangle Triangle;

    void SortAdjacentTriangles(int vIndex,
        std::set<std::shared_ptr<Triangle>> const& tAdj,
        std::vector<std::shared_ptr<Triangle>>& tAdjSorted);

    void CreateSphericalArcs(VETManifoldMesh& mesh, std::multiset<SphericalArc>& arcs);
    void CreateSphericalBisectors(VETManifoldMesh& mesh, std::multiset<SphericalArc>& arcs);
    void CreateBSPTree(std::multiset<SphericalArc>& arcs);
    void InsertArc(SphericalArc const& arcs);

    // Lookup table for indexing into mFaceNormals.
    std::map<std::shared_ptr<Triangle>, int> mTriToNormal;

    // Fixed-size storage for the BSP nodes.
    std::vector<SphericalArc> mNodes;
    int mTreeDepth;
};

template <typename Real>
ExtremalQuery3BSP<Real>::ExtremalQuery3BSP(Polyhedron3<Real> const& polytope)
    :
    ExtremalQuery3<Real>(polytope)
{
    // Create the adjacency information for the polytope.
    VETManifoldMesh mesh;
    auto const& indices = this->mPolytope.GetIndices();
    int const numTriangles = static_cast<int>(indices.size() / 3);
    for (int t = 0; t < numTriangles; ++t)
    {
        int V[3];
        for (int j = 0; j < 3; ++j)
        {
            V[j] = indices[3 * t + j];
        }
        auto triangle = mesh.Insert(V[0], V[1], V[2]);
        mTriToNormal.insert(std::make_pair(triangle, t));
    }

    // Create the set of unique arcs which are used to create the BSP tree.
    std::multiset<SphericalArc> arcs;
    CreateSphericalArcs(mesh, arcs);

    // Create the BSP tree to be used in the extremal query.
    CreateBSPTree(arcs);
}

template <typename Real>
void ExtremalQuery3BSP<Real>::GetExtremeVertices(Vector3<Real> const& direction,
    int& positiveDirection, int& negativeDirection)
{
    // Do a nonrecursive depth-first search of the BSP tree to determine
    // spherical polygon contains the incoming direction D.  Index 0 is the
    // root of the BSP tree.
    int current = 0;
    while (current >= 0)
    {
        SphericalArc& node = mNodes[current];
        int sign = Function<Real>::ISign(Dot(direction, node.normal));
        if (sign >= 0)
        {
            current = node.posChild;
            if (current == -1)
            {
                // At a leaf node.
                positiveDirection = node.posVertex;
            }
        }
        else
        {
            current = node.negChild;
            if (current == -1)
            {
                // At a leaf node.
                positiveDirection = node.negVertex;
            }
        }
    }

    // Do a nonrecursive depth-first search of the BSP tree to determine
    // spherical polygon contains the reverse incoming direction -D.
    current = 0;  // the root of the BSP tree
    while (current >= 0)
    {
        SphericalArc& node = mNodes[current];
        int sign = Function<Real>::ISign(Dot(direction, node.normal));
        if (sign <= 0)
        {
            current = node.posChild;
            if (current == -1)
            {
                // At a leaf node.
                negativeDirection = node.posVertex;
            }
        }
        else
        {
            current = node.negChild;
            if (current == -1)
            {
                // At a leaf node.
                negativeDirection = node.negVertex;
            }
        }
    }
}

template <typename Real> inline
int ExtremalQuery3BSP<Real>::GetNumNodes() const
{
    return static_cast<int>(mNodes.size());
}

template <typename Real> inline
int ExtremalQuery3BSP<Real>::GetTreeDepth() const
{
    return mTreeDepth;
}

template <typename Real>
void ExtremalQuery3BSP<Real>::SortAdjacentTriangles(int vIndex,
    std::set<std::shared_ptr<Triangle>> const& tAdj,
    std::vector<std::shared_ptr<Triangle>>& tAdjSorted)
{
    // Copy the set of adjacent triangles into a vector container.
    int const numTriangles = static_cast<int>(tAdj.size());
    tAdjSorted.resize(tAdj.size());

    // Traverse the triangles adjacent to vertex V using edge-triangle adjacency
    // information to produce a sorted array of adjacent triangles.
    auto tri = *tAdj.begin();
    for (int i = 0; i < numTriangles; ++i)
    {
        for (int prev = 2, curr = 0; curr < 3; prev = curr++)
        {
            if (tri->V[curr] == vIndex)
            {
                tAdjSorted[i] = tri;
                tri = tri->T[prev].lock();
                break;
            }
        }
    }
}

template <typename Real>
void ExtremalQuery3BSP<Real>::CreateSphericalArcs(VETManifoldMesh& mesh,
    std::multiset<SphericalArc>& arcs)
{
    int const prev[3] = { 2, 0, 1 };
    int const next[3] = { 1, 2, 0 };

    for (auto const& element : mesh.GetEdges())
    {
        auto edge = element.second;

        SphericalArc arc;
        arc.nIndex[0] = mTriToNormal[edge->T[0].lock()];
        arc.nIndex[1] = mTriToNormal[edge->T[1].lock()];
        arc.separation = 1;
        arc.normal = Cross(this->mFaceNormals[arc.nIndex[0]], this->mFaceNormals[arc.nIndex[1]]);

        auto adj = edge->T[0].lock();
        int j;
        for (j = 0; j < 3; ++j)
        {
            if (adj->V[j] != edge->V[0] && adj->V[j] != edge->V[1])
            {
                arc.posVertex = adj->V[prev[j]];
                arc.negVertex = adj->V[next[j]];
                break;
            }
        }
        LogAssert(j < 3, "Unexpected condition.");

        arcs.insert(arc);
    }

    CreateSphericalBisectors(mesh, arcs);
}

template <typename Real>
void ExtremalQuery3BSP<Real>::CreateSphericalBisectors(VETManifoldMesh& mesh,
    std::multiset<SphericalArc>& arcs)
{
    std::queue<std::pair<int,int>> queue;
    for (auto const& element : mesh.GetVertices())
    {
        // Sort the normals into a counterclockwise spherical polygon when
        // viewed from outside the sphere.
        auto vertex = element.second;
        int const vIndex = vertex->V;
        std::vector<std::shared_ptr<Triangle>> tAdjSorted;
        SortAdjacentTriangles(vIndex, vertex->TAdjacent, tAdjSorted);
        int const numTriangles = static_cast<int>(vertex->TAdjacent.size());
        queue.push(std::make_pair(0, numTriangles));
        while (!queue.empty())
        {
            std::pair<int, int> item = queue.front();
            queue.pop();
            int i0 = item.first, i1 = item.second;
            int separation = i1 - i0;
            if (separation > 1 && separation != numTriangles - 1)
            {
                if (i1 < numTriangles)
                {
                    SphericalArc arc;
                    arc.nIndex[0] = mTriToNormal[tAdjSorted[i0]];
                    arc.nIndex[1] = mTriToNormal[tAdjSorted[i1]];
                    arc.separation = separation;

                    arc.normal = Cross(this->mFaceNormals[arc.nIndex[0]],
                        this->mFaceNormals[arc.nIndex[1]]);

                    arc.posVertex = vIndex;
                    arc.negVertex = vIndex;
                    arcs.insert(arc);
                }
                int imid = (i0 + i1 + 1) / 2;
                if (imid != i1)
                {
                    queue.push(std::make_pair(i0, imid));
                    queue.push(std::make_pair(imid, i1));
                }
            }
        }
    }
}

template <typename Real>
void ExtremalQuery3BSP<Real>::CreateBSPTree(std::multiset<SphericalArc>& arcs)
{
    // The tree has at least a root.
    mTreeDepth = 1;

    for (auto const& arc : gte::reverse(arcs))
    {
        InsertArc(arc);
    }

    // The leaf nodes are not counted in the traversal of InsertArc.  The
    // depth must be incremented to account for leaves.
    ++mTreeDepth;
}

template <typename Real>
void ExtremalQuery3BSP<Real>::InsertArc(SphericalArc const& arc)
{
    // The incoming arc is stored at the end of the nodes array.
    if (mNodes.size() > 0)
    {
        // Do a nonrecursive depth-first search of the current BSP tree to
        // place the incoming arc.  Index 0 is the root of the BSP tree.
        std::stack<int> candidates;
        candidates.push(0);
        while (!candidates.empty())
        {
            int current = candidates.top();
            candidates.pop();
            SphericalArc* node = &mNodes[current];

            int sign0;
            if (arc.nIndex[0] == node->nIndex[0] || arc.nIndex[0] == node->nIndex[1])
            {
                sign0 = 0;
            }
            else
            {
                Real dot = Dot(this->mFaceNormals[arc.nIndex[0]], node->normal);
                sign0 = Function<Real>::ISign(dot);
            }

            int sign1;
            if (arc.nIndex[1] == node->nIndex[0] || arc.nIndex[1] == node->nIndex[1])
            {
                sign1 = 0;
            }
            else
            {
                Real dot = Dot(this->mFaceNormals[arc.nIndex[1]], node->normal);
                sign1 = Function<Real>::ISign(dot);
            }

            int doTest = 0;
            if (sign0 * sign1 < 0)
            {
                // The new arc straddles the current arc, so propagate it
                // to both child nodes.
                doTest = 3;
            }
            else if (sign0 > 0 || sign1 > 0)
            {
                // The new arc is on the positive side of the current arc.
                doTest = 1;
            }
            else if (sign0 < 0 || sign1 < 0)
            {
                // The new arc is on the negative side of the current arc.
                doTest = 2;
            }
            // else: sign0 = sign1 = 0, in which case no propagation is
            // needed because the current BSP node will handle the correct
            // partitioning of the arcs during extremal queries.

            int depth;

            if (doTest & 1)
            {
                if (node->posChild != -1)
                {
                    candidates.push(node->posChild);
                    depth = static_cast<int>(candidates.size());
                    if (depth > mTreeDepth)
                    {
                        mTreeDepth = depth;
                    }
                }
                else
                {
                    node->posChild = static_cast<int>(mNodes.size());
                    mNodes.push_back(arc);

                    // The push_back can cause a reallocation, so the current
                    // pointer must be refreshed.
                    node = &mNodes[current];
                }
            }

            if (doTest & 2)
            {
                if (node->negChild != -1)
                {
                    candidates.push(node->negChild);
                    depth = static_cast<int>(candidates.size());
                    if (depth > mTreeDepth)
                    {
                        mTreeDepth = depth;
                    }
                }
                else
                {
                    node->negChild = static_cast<int>(mNodes.size());
                    mNodes.push_back(arc);
                }
            }
        }
    }
    else
    {
        // root node
        mNodes.push_back(arc);
    }
}

template <typename Real>
ExtremalQuery3BSP<Real>::SphericalArc::SphericalArc()
    :
    posChild(-1),
    negChild(-1)
{
}

template <typename Real>
bool ExtremalQuery3BSP<Real>::SphericalArc::operator<(SphericalArc const& arc) const
{
    return separation < arc.separation;
}

}
