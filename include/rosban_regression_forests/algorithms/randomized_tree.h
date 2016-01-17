#pragma once

#include "rosban_regression_forests/approximations/approximation_type.h"
#include "rosban_regression_forests/core/training_set.h"
#include "rosban_regression_forests/core/regression_tree.h"
#include "rosban_regression_forests/core/regression_forest.h"

/**
 * Based on 'Extremely Randomized Trees' (Geurts06)
 */
namespace regression_forests
{
namespace RandomizedTree
{
class ExtraTreesConfig
{
public:
  size_t k;
  size_t nMin;
  size_t nbTrees;
  double minVar;
  bool bootstrap;
  ApproximationType apprType;

  ExtraTreesConfig();
  std::vector<std::string> names() const;
  std::vector<std::string> values() const;
  void load(const std::vector<std::string> &names, const std::vector<std::string> &values);
};

double evalSplitScore(const TrainingSet &ls, const TrainingSet::Subset &samples, const OrthogonalSplit &split,
                      enum ApproximationType apprType);
/**
 * k: the number of dimensions used for randomCut
 * nmin: the minimal number of samples per leaf
 * minVariance: if variance is lower than the given threshold,
 *              do not split any further
 */
std::unique_ptr<RegressionTree> learn(const TrainingSet &ls, size_t k, size_t nmin, double minVariance = 0,
                                      enum ApproximationType apprType = ApproximationType::PWC);

std::unique_ptr<RegressionForest> extraTrees(const TrainingSet &ls, size_t k, size_t nmin, size_t nbTrees,
                                             double minVariance = 0, bool bootstrap = false,
                                             enum ApproximationType apprType = ApproximationType::PWC);
};
}
}