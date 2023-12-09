/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <set>

/* HEADER */
#include "ConvexHull.h"

ConvexHull::Facet::Facet(std::vector<size_t>&& _vertexIndices, const bool _isNewFacet) :
  vertexIndices(std::move(_vertexIndices)),
  isNewFacet(_isNewFacet)
{
  std::ranges::sort(vertexIndices);
  neighbors.reserve(vertexIndices.size());
}


ConvexHull::ConvexHull(const Points& unperturbedPoints, const double perturbation, const bool printStats) :
  dimension_(unperturbedPoints.empty() ? 0 : unperturbedPoints.front().size()),
  printStats_(printStats)
{
  powersOfTwelve_.resize(dimension_);
  for (size_t i = 0; i < dimension_; ++i) {
    powersOfTwelve_.at(dimension_ - i - 1) = std::pow(i + 1, 12);
  }
  computeConvexHull_(unperturbedPoints, perturbation);
}

const std::vector<std::vector<size_t>> ConvexHull::getVertexIndices() const
{
  return vertexIndices_;
}

void ConvexHull::computeConvexHull_(const Points& unperturbedPoints, const double perturbation)
{
  // Set seed to ensure deterministic behavior,
  // but use different seeds for different perturbations
  throwExceptionIfInvalidPerturbation_(perturbation, unperturbedPoints);
  const auto* perturbationPtr = &perturbation;
  const auto seed = *reinterpret_cast<const unsigned int*>(perturbationPtr);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> randomDistribution(0.0, 1.0);

  // Perform perturbation of the input points
  points_.clear();
  uniquePointIndexPairs_.resize(unperturbedPoints.size());
  {
    size_t index = 0;
    std::ranges::transform(unperturbedPoints, uniquePointIndexPairs_.begin(), [&index] (const auto& p) { return std::make_pair(p, index++); });
    std::ranges::sort(uniquePointIndexPairs_);
    uniquePointIndexPairs_.erase(std::ranges::unique(uniquePointIndexPairs_, [] (const auto& p1, const auto& p2) { return p1.first == p2.first; }).begin(), uniquePointIndexPairs_.end());
    points_.resize(uniquePointIndexPairs_.size());
    std::ranges::transform(uniquePointIndexPairs_, points_.begin(), [] (const auto& p) { return p.first; });
    for (auto& v : points_) {
      std::ranges::transform(v, v.begin(), [&] (const auto d) { return d + perturbation * randomDistribution(gen) ; });
    }
  }
  throwExceptionIfTooFewPoints_();
  throwExceptionIfNotAllPointsHaveCorrectDimension_();

  // Create simplex using the extreme points
  setInitialSimplex_();

  try {
    // Compute the convex hull for the set of all points using the seed polytope
    growConvexHull_();
  }
  catch (std::invalid_argument& e) {
    // Failed to grow the convex hull. Change perturbation and retry
    std::cerr << e.what() << std::endl;
    const auto newPerturbation = perturbation == 0.0 ? 1e-9 : 10 * perturbation;
    return computeConvexHull_(unperturbedPoints, newPerturbation);
  }

  // Assert that the facets' sets of vertices have the correct dimension
  throwExceptionIfNotAllFacetsFullDimensional_();

  // Construct vector of vertex indices for the facets of the convex hull
  vertexIndices_ = std::vector<std::vector<size_t>>(facets_.size());
  size_t fi = 0;
  for (auto& facet : facets_) {
    vertexIndices_[fi].swap(facet.vertexIndices);
    std::ranges::transform(vertexIndices_[fi], vertexIndices_[fi].begin(), [&] (const auto vi) { return uniquePointIndexPairs_[vi].second; });
    ++fi;
  }
  points_.clear();
  uniquePointIndexPairs_.clear();
  if (printStats_) {
    std::cerr << "Number of distance tests: " << distanceTests << std::endl;
    std::cerr << "Number of hyperplanes created: " << hyperPlanes << std::endl;
  }
}

void ConvexHull::growConvexHull_()
{
  // Check that the input data is correct
  throwExceptionIfNotAllFacetsFullDimensional_();
  throwExceptionIfFacetsUseNonExistingVertices_();

  // Compute origin as the mean of the center points of the seed facets
  updateFacetCenterPoints_();
  const auto origin = computeOrigin_();

  // Compute inwards-oriented facet normals
  updateFacetNormalAndOffset_(origin, facets_);
  throwExceptionIfNotConvexPolytope_();

  // Assign each outer point to a facet
  initializeOutsideSets_();

  // Create a list of all facets that have outside points
  std::vector<FacetIt> facetsWithOutsidePoints;
  facetsWithOutsidePoints.reserve(facets_.size());
  for (auto fIt = facets_.begin(); fIt != facets_.end(); ++fIt) {
    if (!fIt->outsideIndices.empty()) {
      facetsWithOutsidePoints.push_back(fIt);
    }
  }

  // Create new facets using the outside points
  std::ranges::sort(facetsWithOutsidePoints, [] (const auto& a, const auto& b) { return a->farthestOutsidePointDistance < b->farthestOutsidePointDistance; });
  std::vector<FacetIt> visibleFacets;
  std::vector<FacetIt> newFacets;
  while (!facetsWithOutsidePoints.empty()) {
    if (facetsWithOutsidePoints.back()->outsideIndices.empty() ||
        facetsWithOutsidePoints.back()->visible) {
      facetsWithOutsidePoints.pop_back();
      continue;
    }
    auto facetIt = facetsWithOutsidePoints.back();

    // From the outside set of the current facet, find the farthest point
    const auto apexIndex = getAndEraseFarthestPointFromOutsideSet_(*facetIt);

    // Find the set of facets that are visible from the point to be added
    const auto newVisibleFacetsStartIndex = visibleFacets.size();
    std::vector<std::pair<FacetIt, FacetIt>> horizon; // visible-invisible neighboring facet pairs
    getVisibleFacets_(points_.at(apexIndex), facetIt, visibleFacets, horizon);

    // Get the outside points from the visible facets
    std::vector<std::vector<size_t>> unassignedPointIndices(visibleFacets.size() - newVisibleFacetsStartIndex);
    for (size_t i = 0; i < unassignedPointIndices.size(); ++i) {
      unassignedPointIndices[i].swap(visibleFacets[newVisibleFacetsStartIndex + i]->outsideIndices);
    }

    // Create new facets from the apex
    newFacets.clear();
    prepareNewFacets_(apexIndex, horizon, visibleFacets, newFacets);
    connectNeighbors_(apexIndex, horizon, newFacets);
    assert(std::ranges::all_of(newFacets, [dimension = dimension_] (const auto& newFacet) { return newFacet->neighbors.size() == dimension; }));

    // Update the facet normals and offsets
    updateFacetNormalAndOffset_(origin, newFacets);

    // Assign the points belonging to visible facets to the newly created facets
    updateOutsideSets_(unassignedPointIndices, newFacets);

    // Add the new facets with outside points to the vector of all facets with outside points
    for (auto& newFacet : newFacets) {
      newFacet->isNewFacet = false;
      newFacet->visitIndex = (size_t)-1;
    }
    std::erase_if(newFacets, [] (const FacetIt& f) { return f->outsideIndices.empty(); });

    std::ranges::sort(newFacets, [] (const auto& a, const auto& b) { return a->farthestOutsidePointDistance < b->farthestOutsidePointDistance; });
    if (newFacets.empty() ||
        facetsWithOutsidePoints.empty() ||
        newFacets.back()->farthestOutsidePointDistance > facetsWithOutsidePoints.back()->farthestOutsidePointDistance) {
      facetsWithOutsidePoints.insert(facetsWithOutsidePoints.end(), newFacets.begin(), newFacets.end());
    }
    else {
      // facetsWithOutsidePoints.back() has farther point than newFacetsWithOutsidePoints.back()
      newFacets.insert(newFacets.end(), facetsWithOutsidePoints.begin(), facetsWithOutsidePoints.end());
      facetsWithOutsidePoints.swap(newFacets);
    }
  }
  for (const auto& v : visibleFacets) {
    facets_.erase(v);
  }

  // Assert that the facets constitute a convex polytope
  updateFacetCenterPoints_();
  throwExceptionIfNotConvexPolytope_();
}

std::vector<size_t> ConvexHull::getInitialPointIndices_(const Points& points)
{
  const auto dimension = points.size() > 0 ? points[0].size() : 0;
  std::set<size_t> startPointIndexSet;
  for (size_t d = 0; d < dimension; ++d) {
    const auto [min, max] = std::ranges::minmax_element(points, [d] (const auto& a, const auto& b) { return a[d] < b[d]; });
    startPointIndexSet.insert(std::distance(points.cbegin(), min));
    startPointIndexSet.insert(std::distance(points.cbegin(), max));
  }

  for (size_t i = 0; startPointIndexSet.size() <= dimension && i < points.size(); ++i) {
    startPointIndexSet.insert(i);
  }
  return {startPointIndexSet.cbegin(), startPointIndexSet.cend()};
}

std::vector<std::tuple<double, size_t, size_t>> ConvexHull::getPairwiseSquaredDistances_(const std::vector<size_t>& indices, const Points& points)
{
  std::vector<std::tuple<double, size_t, size_t>> distances;
  distances.reserve(indices.size() * (indices.size() + 1) / 2);
  for (size_t i = 0; i < indices.size(); ++i) {
    const auto si = indices[i];
    const auto& pointI = points.at(si);
    for (size_t j = i + 1; j < indices.size(); ++j) {
      const auto sj = indices[j];
      const auto& pointJ = points.at(sj);
      const auto distance = std::transform_reduce(
        pointI.cbegin(), pointI.cend(), pointJ.cbegin(), 0.0,
        std::plus{}, [](auto a, auto b) { return (a - b) * (a - b); }
      );
      distances.emplace_back(distance, si, sj);
    }
  }
  std::ranges::sort(distances);
  return distances;
}

void ConvexHull::setInitialSimplex_()
{
  facets_.clear();
  const auto startPointIndices = getInitialPointIndices_(points_);
  const auto distances = getPairwiseSquaredDistances_(startPointIndices, points_);

  std::set<size_t> indexSet;
  for (size_t di = distances.size(); indexSet.size() <= dimension_; --di) {
    const auto [_, i, j] = distances.at(di - 1);
    indexSet.insert(i);
    indexSet.insert(j);
  }

  // Create initial simplex using the (dimension + 1) first points.
  // The facets have vertices [0, ..., dimension - 1], [1, ..., dimension], ..., [dimension, 0, ..., dimension - 2] in sortedIndices.
  std::vector<size_t> indices(indexSet.cbegin(), indexSet.cend());
  for (size_t i = 0; i <= dimension_; ++i) {
    std::vector<size_t> vertexIndices(dimension_);
    for (size_t j = 0; j < dimension_; ++j) {
      vertexIndices[j] = indices[(i + j) % (dimension_ + 1)];
    }
    facets_.emplace_back(std::move(vertexIndices), false);
  }

  // Update the facets' neighbors
  for (auto fIt1 = facets_.begin(); fIt1 != facets_.end(); ++fIt1) {
    for (auto fIt2 = facets_.begin(); fIt2 != facets_.end(); ++fIt2) {
      if (fIt1 != fIt2) {
        fIt1->neighbors.push_back(fIt2);
      }
    }
  }
}

void ConvexHull::updateFacetCenterPoints_()
{
  const auto oneOverDimension = 1.0 / dimension_;
  for (auto& f : facets_) {
    auto& center = f.center;
    if (center.size() == dimension_) {
      continue;
    }
    center.assign(dimension_, 0.0);
    for (const auto& vi : f.vertexIndices) {
      std::ranges::transform(center, points_.at(vi), center.begin(), std::plus{});
    }
    std::ranges::transform(center, center.begin(), [oneOverDimension] (const auto& c) { return c * oneOverDimension; });
  }
}

std::vector<double> ConvexHull::computeOrigin_() const
{
  assert(!facets_.empty());
  std::vector<double> origin(dimension_, 0.0);
  for (const auto& f : facets_) {
    std::ranges::transform(origin, f.center, origin.begin(), std::plus{});
  }
  std::ranges::transform(origin, origin.begin(), [d = facets_.size()] (const auto& x) { return x / d; });
  return origin;
}

template<template<class T, class All = std::allocator<T>> class Container, class U>
void ConvexHull::updateFacetNormalAndOffset_(const Point& origin, Container<U>& facets)
{
  assert(!points_.empty());
  auto& A = preallocatedA_;
  A.resize(dimension_);
  for (auto& a : A) {
    a.resize(dimension_);
  }
  for (auto& fIt : facets) {
    auto& facet = getFacet_(fIt);
    const auto& firstPoint = points_.at(facet.vertexIndices.front());

    for (size_t i = 0; i + 1 < dimension_; ++i) {
      for (size_t j = 0; j < dimension_; ++j) {
        A[i][j] = points_[facet.vertexIndices[i + 1]][j] - firstPoint[j];
      }
    }
    auto& b = facet.normal;
    b.assign(dimension_, 0.0);
    b.back() = 1.0;
    A.back() = b;

    // Solve A x = b
    overwritingSolveLinearSystemOfEquations_(A, b);
    const auto absSum = 1.0 / std::accumulate(b.cbegin(), b.cend(), 0.0, [](const auto s, const auto& bi) { return s + fabs(bi); });
    std::ranges::transform(b, b.begin(), [absSum] (const auto& bi) { return bi * absSum; });

    facet.offset = scalarProduct_(facet.normal, points_.at(facet.vertexIndices.front()));
    hyperPlanes++;

    // Orient normal inwards
    if (isFacetVisibleFromPoint_(facet, origin)) {
      std::ranges::transform(facet.normal, facet.normal.begin(), [] (const auto d) { return -d; });
      facet.offset = -facet.offset;
    }
  }
}

void ConvexHull::initializeOutsideSets_()
{
  std::set<size_t> vertexIndices;
  for (const auto& facet : facets_) {
    vertexIndices.insert(facet.vertexIndices.cbegin(), facet.vertexIndices.cend());
  }
  for (size_t pi = 0; pi < points_.size(); ++pi) {
    if (vertexIndices.find(pi) == vertexIndices.end()) {
      const auto& point = points_[pi];
      FacetIt farthestFacetIt;
      auto maxDistance = 0.0;
      for (FacetIt fIt = facets_.begin(); fIt != facets_.end(); ++fIt) {
        const auto distance = distance_(*fIt, point);
        if (distance > maxDistance) {
          maxDistance = distance;
          farthestFacetIt = fIt;
        }
      }
      if (maxDistance > 0.0) {
        if (maxDistance > farthestFacetIt->farthestOutsidePointDistance) {
          farthestFacetIt->farthestOutsidePointDistance = maxDistance;
          farthestFacetIt->farthestOutsidePointIndex = farthestFacetIt->outsideIndices.size();
        }
        farthestFacetIt->outsideIndices.push_back(pi);
      }
    }
  }
}

size_t ConvexHull::getAndEraseFarthestPointFromOutsideSet_(Facet& facet)
{
  const auto farthestPointIndex = facet.outsideIndices.at(facet.farthestOutsidePointIndex);
  facet.outsideIndices.at(facet.farthestOutsidePointIndex) = facet.outsideIndices.back();
  facet.outsideIndices.pop_back();
  return farthestPointIndex;
}

void ConvexHull::getVisibleFacets_(const Point& apex,
                                   FacetIt facetIt,
                                   std::vector<FacetIt>& visibleFacets,
                                   std::vector<std::pair<FacetIt, FacetIt>>& horizon) const
{
  facetIt->visible = true;
  facetIt->visitIndex = 0;
  const auto startIndex = visibleFacets.size();
  visibleFacets.push_back(facetIt);
  for (size_t vi = startIndex; vi < visibleFacets.size(); ++vi) {
    const auto visibleFacetIt = visibleFacets.at(vi);
    const auto& visibleFacet = *visibleFacetIt;
    for (auto neighborIt : visibleFacet.neighbors) {
      auto& neighbor = *neighborIt;

      if (neighbor.visitIndex != 0) {
        if (isFacetVisibleFromPoint_(neighbor, apex)) {
          visibleFacets.push_back(neighborIt);
          neighbor.visible = true;
        }
      }

      if (!neighbor.visible) {
        horizon.emplace_back(visibleFacetIt, neighborIt);
      }
      neighbor.visitIndex = 0;
    }
  }
}

void ConvexHull::prepareNewFacets_(const size_t apexIndex,
                                   const std::vector<std::pair<FacetIt, FacetIt>>& horizon,
                                   std::vector<FacetIt>& visibleFacets,
                                   std::vector<FacetIt>& newFacets)
{
  std::vector<Facet> tmpNewFacets;
  newFacets.reserve(horizon.size());
  tmpNewFacets.reserve(std::min(horizon.size(), visibleFacets.size()));
  for (size_t hi = 0; hi < horizon.size(); ++hi) {
    const auto [visibleFacetIt, obscuredFacetIt] = horizon[hi];
    assert(visibleFacetIt->visible);
    assert(!obscuredFacetIt->visible);

    // The new facet has the joint vertices of its parent, plus the index of the apex
    assert(std::ranges::find(visibleFacetIt->vertexIndices, apexIndex) == visibleFacetIt->vertexIndices.end());
    assert(std::ranges::find(obscuredFacetIt->vertexIndices, apexIndex) == obscuredFacetIt->vertexIndices.end());
    std::vector<size_t> vertexIndices(dimension_);
    const auto result = std::ranges::set_intersection(visibleFacetIt->vertexIndices, obscuredFacetIt->vertexIndices, vertexIndices.begin());
    assert(size_t(std::distance(vertexIndices.begin(), result.out) + 1) == dimension_);
    vertexIndices.back() = apexIndex;
    if (hi < visibleFacets.size()) {
      tmpNewFacets.emplace_back(std::move(vertexIndices));
      newFacets.push_back(*(visibleFacets.end() - hi - 1));
    }
    else {
      facets_.emplace_back(std::move(vertexIndices));
      newFacets.push_back(std::prev(facets_.end()));
    }
  }

  // Reuse the space of visible facets, which are to be removed
  for (auto& facet : tmpNewFacets) {
    *visibleFacets.back() = facet;
    visibleFacets.pop_back();
  }

  // Connect new facets to their neighbors
  for (size_t hi = 0; hi < horizon.size(); ++hi) {
    auto newFacetIt = newFacets[hi];
    // The new facet is neighbor to its obscured parent, and vice versa
    auto [visibleFacetIt, obscuredFacetIt] = horizon[hi];
    auto fItIt = std::ranges::find(obscuredFacetIt->neighbors, visibleFacetIt);
    assert(fItIt != obscuredFacetIt->neighbors.end());
    *fItIt = newFacetIt;
    obscuredFacetIt->visitIndex = (size_t)-1;
    newFacetIt->neighbors.push_back(obscuredFacetIt);
  }
}

void ConvexHull::connectNeighbors_(const size_t apexIndex,
                                   const std::vector<std::pair< FacetIt, FacetIt>>& horizon,
                                   std::vector<FacetIt>& newFacets)
{
  auto& peaks = preallocatedPeaks_;
  auto& peakHashes = preallocatedPeakHashes_;
  const auto numOfPeaks = horizon.size() * (dimension_ - 1);
  peaks.resize(numOfPeaks, std::vector<size_t>(dimension_ > 1 ? dimension_ - 2 : 0));
  peakHashes.clear();
  peakHashes.reserve(numOfPeaks);

  size_t peakIndex = 0;
  for (auto& newFacetIt : newFacets) {
    const auto& vertexIndices = newFacetIt->vertexIndices;
    for (const auto& i : vertexIndices) {
      if (i != apexIndex) {
        auto& peak = peaks.at(peakIndex);
        size_t pi = 0;
        for (const auto& j : vertexIndices) {
          if (i != j && j != apexIndex) {
            peak.at(pi++) = j;
          }
        }
        // The vertexIndices are already sorted, so no need to sort peak here.
        const auto hashVal = std::inner_product(peak.cbegin(), peak.cend(), powersOfTwelve_.cbegin(), static_cast<size_t>(0));
        peakHashes.push_back({hashVal, &peak, newFacetIt});
        ++peakIndex;
      }
    }
  }
  /*
  // Simpler, but slower in high dimensions
  std::ranges::sort(peakHashes, [] (const auto& a, const auto& b) {
    const auto& [hash1, peak1, facetIt1] = a;
    const auto& [hash2, peak2, facetIt2] = b;
    return hash1 != hash2 ? hash1 < hash2 : *peak1 < *peak2;
  });

*/
  std::ranges::sort(peakHashes, [] (const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
  // If more than two peaks have the same hash values, it is necessary to sort them by the full vectors
  for (size_t i = 0; i < peakHashes.size();) {
    const auto iVal = std::get<0>(peakHashes[i]);
    size_t j = i + 1;
    while (j < peakHashes.size() && iVal == std::get<0>(peakHashes[j])) {
      ++j;
    }
    if (j > i + 2) {
      // More than two peaks have the same hash value
      std::sort(peakHashes.begin() + i, peakHashes.begin() + j,
                [] (const auto& a, const auto& b) { return *std::get<1>(a) < *std::get<1>(b); });
    }
    i = j;
  }

  // Update neighbors
  for (size_t ri = 0; ri + 1 < peakHashes.size(); ri += 2) {
    const auto& [hash1, peak1, facetIt1] = peakHashes[ri];
    const auto& [hash2, peak2, facetIt2] = peakHashes[ri + 1];

    if ( hash1 != hash2 && *peak1 != *peak2 ) {
      throw std::invalid_argument("Peaks must occur in pairs.");
    }
    facetIt1->neighbors.push_back(facetIt2);
    facetIt2->neighbors.push_back(facetIt1);
  }
}

ConvexHull::Facet* ConvexHull::assignPointToFarthestFacet_(Facet* facet,
                                                           double bestDistance,
                                                           const size_t pointIndex,
                                                           const Point& point,
                                                           const size_t visitIndex) const
{
  // Found a facet for which the point is an outside point
  // Recursively check whether its neighbors are even farther
  facet->visitIndex = visitIndex;
  auto checkNeighbors = true;
  while (checkNeighbors) {
    checkNeighbors = false;
    for (auto& facetIt : facet->neighbors) {
      auto& neighbor = *facetIt;
      if (!neighbor.isNewFacet || neighbor.visitIndex == visitIndex) {
        continue;
      }
      neighbor.visitIndex = visitIndex;
      const auto distance = distance_(neighbor, point);
      if (distance > bestDistance) {
        bestDistance = distance;
        facet = &neighbor;
        checkNeighbors = true;
        break;
      }
    }
  }

  if (bestDistance > facet->farthestOutsidePointDistance) {
    facet->farthestOutsidePointDistance = bestDistance;
    facet->farthestOutsidePointIndex = facet->outsideIndices.size();
  }
  facet->outsideIndices.push_back(pointIndex);
  return facet;
}

void ConvexHull::updateOutsideSets_(const std::vector<std::vector<size_t>>& visibleFacetOutsideIndices,
                                    std::vector<FacetIt>& newFacets) const
{
  assert(!newFacets.empty());

  for (const auto& outsideIndices : visibleFacetOutsideIndices) {
    Facet* facetOfPreviousPoint = nullptr;
    for (size_t pi = 0; pi < outsideIndices.size(); ++pi) {
      const auto pointIndex = outsideIndices.at(pi);
      const auto& point = points_.at(pointIndex);

      if (facetOfPreviousPoint != nullptr) {
        const auto bestDistance = distance_(*facetOfPreviousPoint, point);
        if (bestDistance > 0.0) {
          facetOfPreviousPoint = assignPointToFarthestFacet_(facetOfPreviousPoint, bestDistance, pointIndex, point, pi);
          continue;
        }
      }

      // If the point was not outside the predicted facets, we have to search through all facets
      for (auto& newFacetIt : newFacets) {
        auto& newFacet = *newFacetIt;
        if (facetOfPreviousPoint == &newFacet) {
          continue;
        }
        const auto bestDistance = distance_(newFacet, point);
        if (bestDistance > 0.0) {
          facetOfPreviousPoint = assignPointToFarthestFacet_(&newFacet, bestDistance, pointIndex, point, pi);
          break;
        }
      }
    }
  }
}

bool ConvexHull::isFacetVisibleFromPoint_(const Facet& facet, const Point& point) const
{
  // Returns true if the point is contained in the open negative halfspace of the facet
  return scalarProduct_(facet.normal, point) < facet.offset;
}

double ConvexHull::distance_(const Facet& facet, const Point& point) const
{
  return facet.offset - scalarProduct_(facet.normal, point);
}

double ConvexHull::scalarProduct_(const Point& a, const Point& b) const
{
  ++distanceTests;
  return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), 0.0);
}

void ConvexHull::overwritingSolveLinearSystemOfEquations_(std::vector<std::vector<double> >& A, std::vector<double>& b)
{
  const auto n = A.size();
  assert(n > 0);
  for (auto& ai : A) {
    assert(ai.size() == n);
  }
  assert(b.size() == n);

  // Outer product LU with partial pivoting
  // See Algorithm 3.4.1 in Golub and Van Loan - Matrix Computations, 4th Edition
  for (size_t k = 0; k < n; ++k) {
    // Determine mu with k <= mu < n so abs( A( mu, k ) ) = max( A( k:n-1, k ) )
    auto mu = k;
    auto maxValue = fabs(A[mu][k]);
    for (size_t i = k + 1; i < n; ++i) {
      const auto value = fabs(A[i][k]);
      if (value > maxValue) {
        maxValue = value;
        mu = i;
      }
    }
    if (maxValue == 0.0) {
      throw std::invalid_argument("Singular matrix 1.");
    }
    if (k != mu) {
      A[k].swap(A[mu]);
    }
    // Here, it is utilized that L is not needed
    // (if L is needed, first divide A[ i ][ k ] by A[ k ][ k ], then subtract A[ i ][ k ] * A[ k ][ j ] from A[ i ][ j ])
    const auto& ak = A[k];
    const auto invDiag = 1.0 / ak[k];
    for (size_t i = k + 1; i < n; ++i) {
      auto& ai = A[i];
      const auto factor = ai[k] * invDiag;
      for (size_t j = k + 1; j < n; ++j) {
        ai[j] -= factor * ak[j];
      }
    }
  }
  // LU factorization completed
  // No need to solve Ly = Pb, because b = [0,...,0,1]^T, so y == Pb

  // Solve Ux = y by row-oriented back substitution
  // See Algorithm 3.1.2 in Golub and Van Loan
  for (size_t j = n; j > 0; --j) {
    const auto  i = j - 1;
    const auto sum = std::inner_product(A[i].begin() + j, A[i].end(), b.begin() + j, 0.0);

    if (A[i][i] != 0.0) {
      b[i] = (b[i] - sum) / A[i][i];
    }
    else {
      // Matrix is singular
      if (b[i] == sum) {
        // U(i,i) * x(i) == 0.0 and U(i,i) == 0.0 => x(i) == 0.0 is a solution
        b[i] = 0.0;
      }
      else {
        // U(i,i) * x(i) != 0.0 but U(i,i) == 0.0 => no solution
        throw std::invalid_argument("Singular matrix 2.");
      }
    }
  }
  // b now contains the solution x to Ax = b
}

void ConvexHull::throwExceptionIfNotConvexPolytope_() const
{
  for (const auto& f : facets_) {
    if (std::ranges::any_of(f.neighbors, [&] (const auto& neighbor) { return isFacetVisibleFromPoint_(f, neighbor->center); })) {
      throw std::invalid_argument("Not a convex polytope.");
    }
  }
}

void ConvexHull::throwExceptionIfNotAllFacetsFullDimensional_() const
{
  if (std::ranges::any_of(facets_, [dimension = dimension_] (const auto& f) { return f.vertexIndices.size() != dimension; })) {
    throw std::invalid_argument("All facets must be full dimensional.");
  }
}

void ConvexHull::throwExceptionIfFacetsUseNonExistingVertices_() const
{
  for (const auto& f : facets_) {
    if (std::ranges::any_of(f.vertexIndices, [&] (const auto& vi) { return points_.size() <= vi; } )) {
      throw std::invalid_argument("All facets must consist of existing vertices.");
    }
  }
}

void ConvexHull::throwExceptionIfNotAllPointsHaveCorrectDimension_() const
{
  if (std::ranges::any_of(points_, [dimension = dimension_] (const auto& p) { return p.size() != dimension; })) {
    throw std::invalid_argument("All points must have the correct dimension.");
  }
}

void ConvexHull::throwExceptionIfTooFewPoints_() const
{
  // At least (dimension + 1) points are needed to create a simplex
  if (points_.size() <= dimension_) {
    throw std::invalid_argument("Too few input points to construct convex hull.");
  }
}

void ConvexHull::throwExceptionIfInvalidPerturbation_(const double perturbation, const Points& unperturbedPoints) const
{
  if ( perturbation < 0.0 ) {
    throw std::invalid_argument("Perturbation must be nonnegative.");
  }
  if ( perturbation == 0.0 ) {
    return;
  }
  const auto initialPointIndices = getInitialPointIndices_(unperturbedPoints);
  const auto distances = getPairwiseSquaredDistances_(initialPointIndices, unperturbedPoints);
  const auto maxDistance = std::sqrt(std::get<0>(distances.back()));
  if (perturbation > 0.01 * maxDistance) {
    throw std::invalid_argument("Perturbation is larger than 1 percent of the maximum distance between points." );
  }
}
