#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <cstddef>
#include <list>
#include <vector>

// Computes the vertex indices of the facets constituting the convex hull of the input points
class ConvexHull {
  public:
  using Point = std::vector<double>;
  using Points = std::vector<Point>;

  struct Facet {
    public:
      Facet(std::vector<size_t>&& vertexIndices);
      Facet(const std::vector<size_t>& vertexIndices);

      // indices of the points that are the vertices of the facet
      std::vector<size_t> vertexIndices;

      // the inwards-oriented normal to the facet
      Point normal;

      // the center of the vertices of the facet
      Point center;

      // iterators pointing to the neighboring facets
      std::vector<std::list<Facet>::iterator> neighbors;

      // indices to a set of points that are outside to the facet
      std::vector<size_t> outsideIndices;

      // the offset of the facet from the origin
      double offset = 0.0;

      // indicates whether the facet has been visited in the search for visible facets or in the assignment of outside points
      size_t visitIndex = (size_t)-1;

      // the index of the outside point that is farthest from the facet
      size_t farthestOutsidePointIndex = (size_t)-1;

      // the distance from the farthest point to the facet
      double farthestOutsidePointDistance = 0.0;

      // indicates whether the facet is visible from the point under consideration
      bool visible = false;

      // indicates whether the facet was just created
      bool isNewFacet = true;
  };

    ConvexHull(const Points& points, double perturbation = 0.0);

    const std::vector<std::vector<size_t>> getVertexIndices() const;

    using FacetIt = std::list<Facet>::iterator;

    mutable size_t distanceTests = 0;
    mutable size_t hyperPlanes = 0;

  private:
    void computeConvexHull_(const Points& unperturbedPoints, const double perturbation);

    void growConvexHull_();

    static std::vector<size_t> getInitialPointIndices_(const Points& points);

    static std::vector<std::tuple<double, size_t, size_t>> getPairwiseSquaredDistances_(const std::vector<size_t>& indices, const Points& points);
 
    void setInitialSimplex_();

    void updateFacetCenterPoints_();

    std::vector<double> computeOrigin_() const;

    Facet& getFacet_(Facet& f) {
      return f;
    }
    Facet& getFacet_(FacetIt& fIt) {
      return *fIt;
    }

    template<template<class T, class All = std::allocator<T>> class Container, class U>
    void updateFacetNormalAndOffset_(const Point& origin, Container<U>& facets);

    void initializeOutsideSets_();

    static size_t getAndEraseFarthestPointFromOutsideSet_(Facet& facet);

    void getVisibleFacets_(const Point& apex,
                           FacetIt facetIt,
                           std::vector<FacetIt>& visibleFacets,
                           std::vector<std::pair<FacetIt, FacetIt>>& horizon) const;

    void prepareNewFacets_(size_t apexIndex,
                           const std::vector<std::pair<FacetIt, FacetIt>>& horizon,
                           std::vector<FacetIt>& visibleFacets,
                           std::vector<FacetIt>& newFacets);

    void connectNeighbors_(size_t apexIndex,
                           const std::vector<std::pair<FacetIt, FacetIt>>& horizon,
                           std::vector<FacetIt>& visibleFacets,
                           std::vector<FacetIt>& newFacets);

    Facet* assignPointToFarthestFacet_(Facet* facet,
                                       double bestDistance,
                                       size_t pointIndex,
                                       const Point& point,
                                       size_t visitIndex) const;

    void updateOutsideSets_(const std::vector< std::vector<size_t> >& visibleFacetOutsideIndices,
                            std::vector<FacetIt>& newFacets) const;

    bool isFacetVisibleFromPoint_(const Facet& facet, const Point& point) const;
    double scalarProduct_(const Point& a, const Point& b) const;
    double distance_(const Facet& facet, const Point& point) const;

    static void overwritingSolveLinearSystemOfEquations_(std::vector<std::vector<double> >& A, std::vector<double>& b);


    void throwExceptionIfNotConvexPolytope_() const;
    void throwExceptionIfNotAllFacetsFullDimensional_() const;
    void throwExceptionIfFacetsUseNonExistingVertices_() const;
    void throwExceptionIfInvalidPerturbation_(double perturbation, const Points& unperturbedPoints) const;
    void throwExceptionIfTooFewPoints_() const;
    void throwExceptionIfNotAllPointsHaveCorrectDimension_() const;

    std::vector<std::vector<double>> preallocatedA_;
    std::vector<std::vector<size_t>> preallocatedPeaks_;
    std::vector<std::tuple<size_t, std::vector<size_t>*, FacetIt>> preallocatedPeakHashes_;
    std::vector<size_t> powersOfTwelve_;

    std::vector<std::vector<size_t>> vertexIndices_;
    Points points_;
    std::vector<std::pair<Point, size_t>> uniquePointIndexPairs_;
    std::list<Facet> facets_;
    size_t dimension_ = 0;
};

#endif
