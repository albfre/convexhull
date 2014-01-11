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

      double offset;

      bool visible;

      bool visited;
  };

  // Returns the vertex indices of the facets constituting the convex hull of the input points
  std::vector< std::vector< size_t > > computeConvexHull( const std::vector< std::vector< double > >& points, double perturbation = 0.0 ) /* throw std::invalid_argument */;

  // Computes the convex hull of a set of points given an initial seed hull
  void growConvexHull( const std::vector< std::vector< double > >& points,
                       std::list< Facet >& facets );
}
