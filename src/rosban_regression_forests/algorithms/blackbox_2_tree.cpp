#include "rosban_regression_forests/algorithms/blackbox_2_tree.h"
#include "rosban_regression_forests/approximations/pwc_approximation.h"
#include "rosban_regression_forests/approximations/pwl_approximation.h"

#include "rosban_regression_forests/tools/random.h"
#include "rosban_regression_forests/tools/statistics.h"

#include <iostream>
#include <set>

namespace regression_forests
{
namespace BB2Tree
{
bool SplitEntry::operator<(const SplitEntry &other) const
{
  if (gain == other.gain)
  {
    return node < other.node;
  }
  return gain < other.gain;
}

BB2TreeConfig::BB2TreeConfig()
  : apprType(ApproximationType::PWC), k(1), minPotGain(0), maxLeafs(0), nMin(1), minDensity(0), nbTrees(1)
{
}

std::vector<std::string> BB2TreeConfig::names() const
{
  std::vector<std::string> result = {"ApprType", "k", "minPotGain", "maxLeafs", "nMin", "minDensity", "nbTrees"};
  return result;
}

std::vector<std::string> BB2TreeConfig::values() const
{
  std::vector<std::string> result;
  result.push_back(to_string(apprType));
  result.push_back(std::to_string(k));
  result.push_back(std::to_string(minPotGain));
  result.push_back(std::to_string(maxLeafs));
  result.push_back(std::to_string(nMin));
  result.push_back(std::to_string(minDensity));
  result.push_back(std::to_string(nbTrees));
  return result;
}

void BB2TreeConfig::load(const std::vector<std::string> &colNames, const std::vector<std::string> &colValues)
{
  auto expectedNames = names();
  if (colNames.size() != expectedNames.size())
  {
    throw std::runtime_error("Failed to load extraTreesConfig, mismatch of vector size");
  }
  for (size_t colNo = 0; colNo < colNames.size(); colNo++)
  {
    auto givenName = colNames[colNo];
    auto expectedName = expectedNames[colNo];
    if (givenName.find(expectedName) == std::string::npos)
    {
      throw std::runtime_error("Given name '" + givenName + "' does not match '" + expectedName + "'");
    }
  }
  apprType = loadApproximationType(colValues[0]);
  k = std::stoi(colValues[1]);
  minPotGain = std::stod(colValues[2]);
  maxLeafs = std::stoi(colValues[3]);
  nMin = std::stoi(colValues[4]);
  minDensity = std::stod(colValues[5]);
  nbTrees = std::stoi(colValues[6]);
}

double size(const Eigen::MatrixXd &space)
{
  if (space.cols() != 2)
  {
    throw std::runtime_error("Expecting 2 columns for a space2");
  }
  double prod = 1;
  for (int row = 0; row < space.rows(); row++)
  {
    double dimSize = space(row, 1) - space(row, 0);
    if (dimSize < 0)
    {
      throw std::runtime_error("Space has a max inferior to min");
    }
    prod *= dimSize;
  }
  return prod;
}

// TODO multi-thread version
void populate(TrainingSet &ts, TrainingSet::Subset &samples, const Eigen::MatrixXd &space, int minSize,
              double minDensity, EvalFunc eval)
{
  int minSamplesByDensity = ceil(minDensity * size(space));
  int wishedSamples = std::max(minSize, minSamplesByDensity) - samples.size();
  if (wishedSamples <= 0)
  {
    return;
  }
  std::vector<Eigen::VectorXd> inputs = regression_forests::getUniformSamples(space, wishedSamples);
  for (const Eigen::VectorXd &i : inputs)
  {
    samples.push_back(ts.size());
    ts.push(Sample(i, eval(i)));
  }
}

Approximation *getApproximation(const TrainingSet &ts, const TrainingSet::Subset &samples, ApproximationType apprType)
{
  switch (apprType)
  {
    case PWC:
      return new PWCApproximation(Statistics::mean(ts.values(samples)));
    case PWL:
      return new PWLApproximation(ts.inputs(samples), ts.values(samples));
  }
  throw std::runtime_error("Unknown approximation type in getApproximation");
}

double potentialGain(const TrainingSet &ts, const TrainingSet::Subset &samples, Approximation *a,
                     const Eigen::MatrixXd &space)
{
  double squaredSum = 0;
  for (int s : samples)
  {
    double err = ts(s).getOutput() - a->eval(ts(s).getInput());
    squaredSum += err * err;
  }
  double var = squaredSum / samples.size();
  double res = var * size(space);
  return res;
}

double avgSquaredErrors(const TrainingSet &ls, const TrainingSet::Subset &samples, enum ApproximationType apprType)
{
  switch (apprType)
  {
    case PWC:
      return Statistics::variance(ls.values(samples));
    case PWL:
    {
      std::vector<Eigen::VectorXd> inputs = ls.inputs(samples);
      std::vector<double> outputs = ls.values(samples);
      PWLApproximation a(inputs, outputs);
      double sumSquaredError = 0;
      for (size_t i = 0; i < inputs.size(); i++)
      {
        double error = a.eval(inputs[i]) - outputs[i];
        sumSquaredError += error * error;
      }
      return sumSquaredError / inputs.size();
    }
  }
  throw std::runtime_error("Unknown ApprType");
}

double evalSplitScore(const TrainingSet &ls, const TrainingSet::Subset &samples, const OrthogonalSplit &split,
                      enum ApproximationType apprType)
{
  std::vector<int> samplesUpper, samplesLower;
  ls.applySplit(split, samples, samplesLower, samplesUpper);
  double nbSamples = samples.size();
  // This happened once in Randomized Tree, but no idea why
  if (samplesLower.size() == 0 || samplesUpper.size() == 0)
  {
    std::ostringstream oss;
    oss << "One of the sample Set is empty, this should never happen" << std::endl;
    oss << "Split: (" << split.dim << "," << split.val << ")" << std::endl;
    oss << "Samples:" << std::endl;
    std::vector<double> dimValues = ls.inputs(samples, split.dim);
    std::sort(dimValues.begin(), dimValues.end());
    for (double v : dimValues)
    {
      oss << v << std::endl;
    }
    std::cout << oss.str();
    throw std::logic_error("Empty sample set provided to evalSplitScore");
  }
  double varAll = avgSquaredErrors(ls, samples, apprType);
  double varLower = avgSquaredErrors(ls, samplesLower, apprType);
  double varUpper = avgSquaredErrors(ls, samplesUpper, apprType);
  double weightedNewVar = (varLower * samplesLower.size() + varUpper * samplesUpper.size()) / nbSamples;

  if (varAll == 0)
  {
    return 0;
  }
  return (varAll - weightedNewVar) / varAll;
}

SplitEntry getBestSplitEntry(RegressionNode *node, const TrainingSet &ts, TrainingSet::Subset &samples,
                             const Eigen::MatrixXd &space, int k, int nMin, ApproximationType apprType)
{
  auto generator = regression_forests::get_random_engine();
  std::vector<size_t> dimCandidates = regression_forests::getKDistinctFromN(k, ts.getInputDim(), &generator);
  std::vector<OrthogonalSplit> splitCandidates;
  splitCandidates.reserve(k);
  for (int i = 0; i < k; i++)
  {
    size_t dim = dimCandidates[i];
    ts.sortSubset(samples, dim);
    double sValMin = ts(samples[nMin - 1]).getInput(dim);
    double sValMax = ts(samples[samples.size() - nMin]).getInput(dim);
    if (sValMin == sValMax)
    {
      std::cout << "Similar Values for sValMin and sValMax: " << sValMin << std::endl;
      std::cout << "Space: " << std::endl
                << space << std::endl;
      std::cout << "SpaceSize: " << std::endl
                << size(space) << std::endl;
    }
    std::uniform_real_distribution<double> distribution(sValMin, sValMax);
    splitCandidates.push_back(OrthogonalSplit(dim, distribution(generator)));
  }
  // If no splits are available: throw an exception (if all inputs are similar
  // concerning the k dimension)
  if (splitCandidates.size() == 0)
  {
    throw std::runtime_error("No possible splits");
  }

  // Find best split candidate
  size_t bestSplitIdx = 0;
  double bestSplitScore = evalSplitScore(ts, samples, splitCandidates[0], apprType);
  for (size_t splitIdx = 1; splitIdx < splitCandidates.size(); splitIdx++)
  {
    double splitScore = evalSplitScore(ts, samples, splitCandidates[splitIdx], apprType);
    if (splitScore > bestSplitScore)
    {
      bestSplitScore = splitScore;
      bestSplitIdx = splitIdx;
    }
  }

  SplitEntry e;
  e.node = node;
  ;
  e.gain = bestSplitScore * size(space);
  e.split = splitCandidates[bestSplitIdx];
  e.samples = samples;
  e.space = space;
  return e;
}

void treat(RegressionNode *node, TrainingSet &ts, TrainingSet::Subset &samples, const Eigen::MatrixXd &space,
           const BB2TreeConfig &c, std::set<SplitEntry> &splitCandidates)
{
  populate(ts, samples, space, 2 * c.nMin, c.minDensity, c.eval);
  node->a = getApproximation(ts, samples, c.apprType);
  double potGain = potentialGain(ts, samples, node->a, space);
  if (potGain >= c.minPotGain)
  {
    try
    {
      splitCandidates.insert(getBestSplitEntry(node, ts, samples, space, c.k, c.nMin, c.apprType));
    }
    catch (const std::runtime_error &exc)
    {
      // There was no available split
      std::cout << "No available splits" << std::endl;
    }
  }
}

std::unique_ptr<RegressionTree> bb2Tree(const BB2TreeConfig &config)
{
  std::unique_ptr<RegressionTree> tree(new RegressionTree);
  std::set<SplitEntry> splitCandidates;
  TrainingSet ts(config.space.rows());
  TrainingSet::Subset samples;
  tree->root = new RegressionNode();
  treat(tree->root, ts, samples, config.space, config, splitCandidates);
  size_t nbLeafs = 1;
  while (splitCandidates.size() > 0 && nbLeafs < config.maxLeafs)
  {
    // Retrieving bestSplit Candidates
    SplitEntry entry = *std::prev(splitCandidates.end());
    splitCandidates.erase(std::prev(splitCandidates.end()));
    RegressionNode *node = entry.node;
    // Splitting Samples
    TrainingSet::Subset lSamples, uSamples;
    ts.applySplit(entry.split, entry.samples, lSamples, uSamples);
    // Modifying node
    delete (node->a);
    node->a = NULL;
    node->s = entry.split;
    node->lowerChild = new RegressionNode();
    node->upperChild = new RegressionNode();
    // Getting spaces
    Eigen::MatrixXd lowerSpace = entry.space;
    lowerSpace(entry.split.dim, 1) = entry.split.val;
    Eigen::MatrixXd upperSpace = entry.space;
    upperSpace(entry.split.dim, 0) = entry.split.val;
    // Treating childs
    try
    {
      treat(node->lowerChild, ts, lSamples, lowerSpace, config, splitCandidates);
      treat(node->upperChild, ts, uSamples, upperSpace, config, splitCandidates);
    }
    catch (const std::logic_error &exc)
    {
      std::cerr << "Error while treating one of the child" << std::endl;
      std::cerr << "\tFatherSpace" << std::endl
                << entry.space << std::endl;
      std::cerr << "\tLowerSpace" << std::endl
                << lowerSpace << std::endl;
      std::cerr << "\tUpperSpace" << std::endl
                << upperSpace << std::endl;
      std::cerr << "SplitCandidate.gain" << entry.gain << std::endl;
    }
    // Updating number of leafs
    nbLeafs++;
  }
  std::cout << "NbLeafs: " << nbLeafs << "/" << config.maxLeafs << std::endl;
  return tree;
}

std::unique_ptr<RegressionForest> bb2Forest(const BB2TreeConfig &config)
{
  std::unique_ptr<RegressionForest> forest(new RegressionForest);
  for (int i = 0; i < config.nbTrees; i++)
  {
    forest->push(bb2Tree(config));
  }
  return forest;
}
}  // Namespace BB2Tree
}
}
