#include "rhoban_regression_forests/core/tree.h"

#include "rhoban_utils/io_tools.h"

#include <limits>
#include <iostream>

namespace regression_forests
{
Tree::Tree() : root(NULL)
{
}

Tree::~Tree()
{
  if (root != NULL)
  {
    delete (root);
  }
}

size_t Tree::nbLeafs() const
{
  return root->nbLeafs();
}

size_t Tree::maxSplitDim() const
{
  return root->maxSplitDim();
}

double Tree::getValue(const Eigen::VectorXd& input) const
{
  return root->getValue(input);
}

Eigen::VectorXd Tree::getGrad(const Eigen::VectorXd& input) const
{
  return root->getGrad(input);
}

double Tree::getMax(const Eigen::MatrixXd& limits) const
{
  return getMaxPair(limits).first;
}

Eigen::VectorXd Tree::getArgMax(const Eigen::MatrixXd& limits) const
{
  return getMaxPair(limits).second;
}

std::pair<double, Eigen::VectorXd> Tree::getMinPair(const Eigen::MatrixXd& limits) const
{
  Eigen::MatrixXd localLimits = limits;
  return root->getMinPair(localLimits);
}

std::pair<double, Eigen::VectorXd> Tree::getMaxPair(const Eigen::MatrixXd& limits) const
{
  Eigen::MatrixXd localLimits = limits;
  return root->getMaxPair(localLimits);
}

void Tree::fillProjection(std::vector<std::vector<Eigen::VectorXd>>& out, Node* currentNode,
                          const std::vector<int>& freeDimensions, Eigen::MatrixXd& limits)
{
  // If leaf, fill vector and return
  if (currentNode->isLeaf())
  {
    out.push_back(currentNode->project(freeDimensions, limits));
    return;
  }
  int splitDim = currentNode->s.dim;
  double splitVal = currentNode->s.val;
  double oldMin = limits(splitDim, 0);
  double oldMax = limits(splitDim, 1);
  // If lowerChild shares a common space with limits
  if (oldMin <= splitVal)
  {
    // set a new limit if necessary
    limits(splitDim, 1) = std::min(splitVal, oldMax);
    // Explore subtree
    fillProjection(out, currentNode->lowerChild, freeDimensions, limits);
    // Restor old limit
    limits(splitDim, 1) = oldMax;
  }
  // If upperChild shares a common space with limits
  if (oldMax > splitVal)
  {
    // set a new limit if necessary
    limits(splitDim, 0) = std::max(splitVal, oldMin);
    // Explore subtree
    fillProjection(out, currentNode->upperChild, freeDimensions, limits);
    // Restor old limit
    limits(splitDim, 0) = oldMin;
  }
}

std::vector<std::vector<Eigen::VectorXd>> Tree::project(const std::vector<int>& freeDimensions,
                                                        const Eigen::MatrixXd& limits)
{
  std::vector<std::vector<Eigen::VectorXd>> out;
  Eigen::MatrixXd localLimits = limits;  // Copying to allow modification
  fillProjection(out, root, freeDimensions, localLimits);
  return out;
}

void Tree::addSubTree(Node* node, Eigen::MatrixXd& limits, const Tree& other, double otherWeight)
{
  if (otherWeight == 0)
  {
    return;
  }
  if (node->isLeaf())
  {
    // Deep copy of the 'other' subTree corresponding to limits
    Node* subTreeRoot = other.root->subTreeCopy(limits);
    // If the subTree is void, then nothing needs to be done
    if (subTreeRoot == NULL)
    {
      return;
    }
    // Add current approximation to all leafs of the subTree
    if (node->a)
    {
      subTreeRoot->addApproximation(node->a, 1.0 / otherWeight);
    }
    // Replace approximation, split and childs by subTreeRoot values
    node->a = subTreeRoot->a;
    node->s = subTreeRoot->s;
    node->lowerChild = subTreeRoot->lowerChild;
    node->upperChild = subTreeRoot->upperChild;
    // Replace the father of subTreeRoot childs if necessary
    if (!node->isLeaf())
    {
      node->lowerChild->father = node;
      node->upperChild->father = node;
    }
    // Set pointers to NULL and delete subTreeRoot
    subTreeRoot->lowerChild = NULL;
    subTreeRoot->upperChild = NULL;
    delete (subTreeRoot);
    return;
  }
  int sDim = node->s.dim;
  double sVal = node->s.val;
  double oldMin = limits(sDim, 0);
  double oldMax = limits(sDim, 1);
  // Add subTree to lowerchild (if lowerchild is in limits)
  if (limits(sDim, 0) <= sVal)
  {
    limits(sDim, 1) = sVal;
    addSubTree(node->lowerChild, limits, other, otherWeight);
    limits(sDim, 1) = oldMax;
  }
  // Add subTree to upperchild (if upperchild is in limits)
  if (limits(sDim, 1) > sVal)
  {
    limits(sDim, 0) = sVal;
    addSubTree(node->upperChild, limits, other, otherWeight);
    limits(sDim, 0) = oldMin;
  }
}

void Tree::avgTree(const Tree& other, double otherWeight)
{
  size_t dim = 1 + std::max(maxSplitDim(), other.maxSplitDim());
  Eigen::MatrixXd limits(dim, 2);
  for (size_t d = 0; d < dim; d++)
  {
    limits(d, 0) = std::numeric_limits<double>::lowest();
    limits(d, 1) = std::numeric_limits<double>::max();
  }
  avgTree(other, otherWeight, limits);
}

void Tree::avgTree(const Tree& other, double otherWeight, const Eigen::MatrixXd& limits)
{
  Eigen::MatrixXd localLimits = limits;
  addSubTree(root, localLimits, other, otherWeight);
}

std::unique_ptr<Tree> Tree::avgTrees(const Tree& t1, const Tree& t2, double w1, double w2,
                                     const Eigen::MatrixXd& limits)
{
  Eigen::MatrixXd localLimits = limits;
  std::unique_ptr<Tree> result(new Tree());
  result->root = new Node();
  Node::parallelMerge(*result->root, *t1.root, *t2.root, w1, w2, localLimits);
  return result;
}

std::unique_ptr<Tree> Tree::project(const Eigen::MatrixXd& limits) const
{
  std::unique_ptr<Tree> t(new Tree);
  t->root = new Node(NULL, NULL);
  Eigen::MatrixXd localLimits = limits;
  t->addSubTree(t->root, localLimits, *this, 1.0);
  return t;
}

void Tree::apply(Eigen::MatrixXd& limits, Node::Function f)
{
  root->apply(limits, f);
}

void Tree::applyOnLeafs(Eigen::MatrixXd& limits, Node::Function f)
{
  root->applyOnLeafs(limits, f);
}

int Tree::write(std::ostream& out) const
{
  char has_root = 1;
  if (root == nullptr)
    has_root = 0;
  int bytes_written = 0;
  bytes_written += rhoban_utils::write<char>(out, has_root);
  if (has_root)
  {
    bytes_written += root->write(out);
  }
  return bytes_written;
}

int Tree::read(std::istream& in)
{
  // First free currently used data
  if (root != nullptr)
    delete (root);
  // Then read
  char has_root;
  int bytes_read = 0;
  bytes_read += rhoban_utils::read<char>(in, &has_root);
  if (has_root)
  {
    root = new Node();
    bytes_read += root->read(in);
  }
  return bytes_read;
}

Tree* Tree::clone() const
{
  Tree* copy = new Tree();
  copy->root = root->clone();
  return copy;
}

}  // namespace regression_forests
