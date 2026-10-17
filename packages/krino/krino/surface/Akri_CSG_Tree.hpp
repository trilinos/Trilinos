// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SURFACE_AKRI_CSG_TREE_HPP_
#define KRINO_KRINO_SURFACE_AKRI_CSG_TREE_HPP_
#include <memory>
#include <string>
#include <vector>

namespace krino {

class Surface;

enum class CSGOperation
{
    Leaf,
    Union,
    Intersection,
    Difference
};

class CSGTreeNode
{
public:
    explicit CSGTreeNode(std::shared_ptr<Surface> surf);

    CSGTreeNode(
        CSGOperation operation,
        std::shared_ptr<CSGTreeNode> lhs,
        std::shared_ptr<CSGTreeNode> rhs);
    CSGTreeNode(
            CSGOperation operation,
            const std::vector<std::shared_ptr<CSGTreeNode>> & operands);

    CSGOperation node_operator() const { return myOp; }
    const std::shared_ptr<Surface> & surface() const;
    const std::shared_ptr<CSGTreeNode> & left_operand() const;
    const std::shared_ptr<CSGTreeNode> & right_operand() const;
    const std::vector<std::shared_ptr<CSGTreeNode>> & vector_operands() const;
private:
    CSGOperation myOp;
    std::shared_ptr<Surface> mySurface; // Used only for Leaf nodes
    std::vector<std::shared_ptr<CSGTreeNode>> myOperands;

};

// ------------------------------------------------------------
// CSGTree: responsible for building and storing the CSG tree
// ------------------------------------------------------------
class CSGTree
{
private:
    std::shared_ptr<CSGTreeNode> myRoot;

    static void gatherSurfacesRecursive(
        const std::shared_ptr<CSGTreeNode>& node,
        std::vector<std::shared_ptr<Surface>>& surfaces);

public:
    CSGTree();

    explicit CSGTree(std::shared_ptr<Surface> surface);

    explicit CSGTree(std::shared_ptr<CSGTreeNode> root);

    static CSGTree unite(const CSGTree& a, const CSGTree& b);
    static CSGTree unite(std::shared_ptr<Surface> a, std::shared_ptr<Surface> b) { return unite(CSGTree(a), CSGTree(b)); }
    static CSGTree unite(const std::vector<CSGTree> & operands);
    static CSGTree unite(const std::vector<std::shared_ptr<Surface>> & surfaces);

    static CSGTree intersect(const CSGTree& a, const CSGTree& b);
    static CSGTree intersect(std::shared_ptr<Surface> a, std::shared_ptr<Surface> b) { return intersect(CSGTree(a), CSGTree(b)); }
    static CSGTree intersect(const std::vector<CSGTree> & operands);
    static CSGTree intersect(const std::vector<std::shared_ptr<Surface>> & surfaces);

    static CSGTree subtract(const CSGTree& a, const CSGTree& b);
    static CSGTree subtract(std::shared_ptr<Surface> a, std::shared_ptr<Surface> b) { return subtract(CSGTree(a), CSGTree(b)); }

    // Mutating builder methods
    void setRoot(std::shared_ptr<Surface> surface);

    void makeUnionWith(const CSGTree& other);
    void makeUnionWith(const std::shared_ptr<Surface> other) { makeUnionWith(CSGTree(other)); }

    void makeIntersectionWith(const CSGTree& other);
    void makeIntersectionWith(const std::shared_ptr<Surface> other) { makeIntersectionWith(CSGTree(other)); }

    void subtractInPlace(const CSGTree& other);
    void subtractInPlace(const std::shared_ptr<Surface> other) { subtractInPlace(CSGTree(other)); }

    bool empty() const;
    void clear();

    std::shared_ptr<CSGTreeNode> root() const;

    std::shared_ptr<CSGTreeNode> requireRoot() const;

    std::vector<std::shared_ptr<Surface>> gather_surfaces() const;

    static std::string operation_name(CSGOperation op);
};


}

#endif /* KRINO_KRINO_SURFACE_AKRI_CSG_TREE_HPP_ */
