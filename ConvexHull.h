#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <cstddef>
#include <list>
#include <vector>

namespace ConvexHull {
  struct Facet {
    public:
      Facet( const std::vector< size_t >& vertexIndices );

      // indices of the points that are the vertices of the facet
      std::vector< size_t > vertexIndices;

      // the inwards-oriented normal to the facet
      std::vector< double > normal;

      // the center of the vertices of the facet
      std::vector< double > center;

      // iterators pointing to the neighboring facets
      std::vector< std::list< Facet >::iterator > neighbors;

      // indices to a set of points that are outside to the facet
      std::vector< size_t > outsideIndices;

      // the offset of the facet from the origin
      double offset;

      // indicates whether the facet has been visited in the search for visible facets or in the assignment of outside points
      size_t visitIndex;

      // the index of the outside point that is farthest from the facet
      size_t farthestOutsidePointIndex;

      // the distance from the farthest point to the facet
      double farthestOutsidePointDistance;

      // indicates whether the facet is visible from the point under consideration
      bool visible;

      // indicates whether the facet was just created
      bool isNewFacet;

      bool hasObscured;
  };

  // Returns the vertex indices of the facets constituting the convex hull of the input points
  std::vector< std::vector< size_t > > computeConvexHull( const std::vector< std::vector< double > >& points, double perturbation = 0.0 ) /* throw std::invalid_argument */;

  // Computes the convex hull of a set of points given an initial seed hull
  void growConvexHull( const std::vector< std::vector< double > >& points,
                       std::list< Facet >& facets );
}

#endif
