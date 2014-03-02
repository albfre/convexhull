/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <set>

/* HEADER */
#include "ConvexHull.h"

// For profiling
#if 0
#define INLINE_ATTRIBUTE __attribute__ ((noinline))
#else
#define INLINE_ATTRIBUTE
#endif


using namespace std;

namespace ConvexHull {

Facet::Facet( const vector< size_t >& _vertexIndices ) :
  vertexIndices( _vertexIndices ),
  visible( false ),
  visitIndex( (size_t)-1 ),
  isNewFacet( true ),
  farthestOutsidePointDistance( 0.0 )
{
  neighbors.reserve( vertexIndices.size() );
}

typedef list< Facet >::iterator FacetIt;
typedef list< Facet >::const_iterator FacetConstIt;

/* Anonymous namespace function declarations */
namespace {
  vector< vector< size_t > > INLINE_ATTRIBUTE computeConvexHull_( const vector< vector< double > >& unperturbedPoints, double perturbation, size_t depth );
  vector< size_t > INLINE_ATTRIBUTE getInitialPolytopeVertexIndices_( const vector< vector< double > >& points );
  void INLINE_ATTRIBUTE getInitialSimplex_( const vector< vector< double > >& points, list< Facet >& facets );
  vector< vector< double > > INLINE_ATTRIBUTE perturbPoints_( const vector< vector< double > >& points, double perturbation );
  void INLINE_ATTRIBUTE updateFacetCenterPoints_( const vector< vector< double > >& points, list< Facet >& facets );
  vector< double > INLINE_ATTRIBUTE computeOrigin_( const list< Facet >& facets );
  template< template< class T, class All = std::allocator< T > > class Container, class U >
  void INLINE_ATTRIBUTE updateFacetNormalAndOffset_( const vector< vector< double > >& points, const vector< double >& origin, Container< U >& facets, vector< vector< double > >& preallocatedA );
  void INLINE_ATTRIBUTE initializeOutsideSets_( const vector< vector< double > >& points, list< Facet >& facets );
  size_t INLINE_ATTRIBUTE getAndEraseFarthestPointFromOutsideSet_( const vector< vector< double > >& points, Facet& facet );
  void INLINE_ATTRIBUTE getVisibleFacets_( const vector< double >& apex, FacetIt facetIt, vector< FacetIt >& visibleFacets, vector< pair< FacetIt, FacetIt > >& horizon );
  size_t INLINE_ATTRIBUTE getHashValue_( const vector< size_t >& v );
  void INLINE_ATTRIBUTE createNewFacets_( size_t apexIndex, const vector< pair< FacetIt, FacetIt > >& horizon, list< Facet >& facets, vector< FacetIt >& visibleFacets, vector< FacetIt >& newFacets, vector< vector< size_t > >& preallocatedPeaks );
  void INLINE_ATTRIBUTE updateOutsideSets_( const vector< vector< double > >& points, const vector< size_t >& unassignedPointIndices, vector< FacetIt >& newFacets );
  vector< size_t > INLINE_ATTRIBUTE getOutsidePointIndicesFromFacets_( const vector< FacetIt >& facets, size_t startIndex );
  bool INLINE_ATTRIBUTE isFacetVisibleFromPoint_( const Facet& facet, const vector< double >& point );
  double INLINE_ATTRIBUTE distance_( const Facet& facet, const vector< double >& point );
  double INLINE_ATTRIBUTE scalarProduct_( const vector< double >& a, const vector< double >& b );
  void INLINE_ATTRIBUTE overwritingSolveLinearSystemOfEquations_( vector< vector< double > >& A, vector< double >& b );

  void INLINE_ATTRIBUTE throwExceptionIfNotConvexPolytope_( const list< Facet >& facets );
  void INLINE_ATTRIBUTE throwExceptionIfNotAllFacetsFullDimensional_( const list< Facet >& facets, size_t dimension );
  void INLINE_ATTRIBUTE throwExceptionIfFacetsUseNonExistingVertices_( const list< Facet >& facets, const vector< vector< double > >& points );
  void INLINE_ATTRIBUTE throwExceptionIfNotAllPointsHaveCorrectDimension_( const vector< vector< double > >& points, size_t dimension );
  void INLINE_ATTRIBUTE throwExceptionIfTooFewPoints_( const vector< vector< double > >& points );
  void INLINE_ATTRIBUTE throwExceptionIfInvalidPerturbation_( double perturbation, const vector< vector< double > >& points );

  struct FirstComparator {
    template< class T, class U >
    bool operator() ( const pair< T, U >& p1,
                      const pair< T, U >& p2 ) const { return p1.first < p2.first; } };

  struct SecondSecondPtrComparator {
    template< class T, class U, class V >
    bool operator() ( const pair< T, pair< U, V* > >& p1,
                      const pair< T, pair< U, V* > >& p2 ) const { return *p1.second.second < *p2.second.second; } };

  struct FarthestPointDistanceComparator {
    bool operator() ( const FacetIt& fIt1,
                      const FacetIt& fIt2 ) const { return fIt1->farthestOutsidePointDistance < fIt2->farthestOutsidePointDistance; } };

  struct DimensionComparator {
    DimensionComparator( size_t dimension ) : dimension_( dimension ) {}
    bool operator() ( const vector< double >& p1, const vector< double >& p2 ) const { return p1[ dimension_ ] < p2[ dimension_ ]; }
    private:
      size_t dimension_;
  };

  struct IsVisiblePredicate { bool operator() ( const FacetIt& f ) const { return f->visible; } };

  template< class T > struct FacetExtractor;

  template<>
  struct FacetExtractor< FacetIt > {
    static Facet& getFacet( vector< FacetIt >::iterator fItIt ) { return *( *fItIt ); }
  };

  template<>
  struct FacetExtractor< Facet > {
    static Facet& getFacet( FacetIt fIt ) { return *fIt; }
  };
}

vector< vector< size_t > > computeConvexHull( const vector< vector< double > >& unperturbedPoints,
                                              double perturbation )
{
  // Always use the same seed to ensure deterministic behavior
  srand( 4711 );

  // Check that the input data is correct
  throwExceptionIfTooFewPoints_( unperturbedPoints );
  const size_t dimension = unperturbedPoints.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( unperturbedPoints, dimension );

  return computeConvexHull_( unperturbedPoints, perturbation, 0 );
}

void growConvexHull( const vector< vector< double > >& points,
                     list< Facet >& facets )
{
  // Check that the input data is correct
  throwExceptionIfTooFewPoints_( points );
  const size_t dimension = points.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( points, dimension );
  throwExceptionIfNotAllFacetsFullDimensional_( facets, dimension );
  throwExceptionIfFacetsUseNonExistingVertices_( facets, points );

  // The vertex indices are assumed to be sorted when new facets are created
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    sort( fIt->vertexIndices.begin(), fIt->vertexIndices.end() );
  }

  // Update the facet center points to the mean of the vertex points
  updateFacetCenterPoints_( points, facets );

  // Compute origin as the mean of the center points of the seed facets
  const vector< double > origin = computeOrigin_( facets );

  // Compute inwards-oriented facet normals
  vector< vector< double > > preallocatedA;
  updateFacetNormalAndOffset_( points, origin, facets, preallocatedA );
  throwExceptionIfNotConvexPolytope_( facets );

  // Assign each outer point to a facet
  initializeOutsideSets_( points, facets );

  // Create a list of all facets that have outside points
  vector< FacetIt > facetsWithOutsidePoints;
  facetsWithOutsidePoints.reserve( facets.size() );
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    if ( fIt->outsideIndices.size() > 0 ) {
      facetsWithOutsidePoints.push_back( fIt );
    }
  }

  // Create new facets using the outside points
  sort( facetsWithOutsidePoints.begin(), facetsWithOutsidePoints.end(), FarthestPointDistanceComparator() );
  vector< FacetIt > visibleFacets;
  vector< vector< size_t > > preallocatedPeaks;
  while ( facetsWithOutsidePoints.size() > 0 ) {
    if ( facetsWithOutsidePoints.back()->outsideIndices.size() == 0 ||
         facetsWithOutsidePoints.back()->visible ) {
      facetsWithOutsidePoints.pop_back();
      continue;
    }
    FacetIt facetIt = facetsWithOutsidePoints.back();

    // From the outside set of the current facet, find the farthest point
    const size_t apexIndex = getAndEraseFarthestPointFromOutsideSet_( points, *facetIt );
    const vector< double >& apex = points[ apexIndex ];

    // Find the set of facets that are visible from the point to be added
    vector< pair< FacetIt, FacetIt > > horizon; // visible-invisible neighboring facet pairs
    size_t newVisibleFacetsStartIndex = visibleFacets.size();
    getVisibleFacets_( apex, facetIt, visibleFacets, horizon );

    // Get the outside points from the visible facets
    const vector< size_t >& unassignedPointIndices = getOutsidePointIndicesFromFacets_( visibleFacets, newVisibleFacetsStartIndex );

    // Create new facets from the apex
    vector< FacetIt > newFacets;
    createNewFacets_( apexIndex, horizon, facets, visibleFacets, newFacets, preallocatedPeaks );

    // Update the facet normals and offsets
    updateFacetNormalAndOffset_( points, origin, newFacets, preallocatedA );

    // Assign the points belonging to visible facets to the newly created facets
    updateOutsideSets_( points, unassignedPointIndices, newFacets );

    // Add the new facets with outside points to the vector of all facets with outside points
    vector< FacetIt > newFacetsWithOutsidePoints;
    newFacetsWithOutsidePoints.reserve( newFacets.size() );
    for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
      newFacets[ ni ]->isNewFacet = false;
      newFacets[ ni ]->visitIndex = (size_t)-1;
      if ( newFacets[ ni ]->outsideIndices.size() > 0 ) {
        newFacetsWithOutsidePoints.push_back( newFacets[ ni ] );
      }
    }
    sort( newFacetsWithOutsidePoints.begin(), newFacetsWithOutsidePoints.end(), FarthestPointDistanceComparator() );
    if ( newFacetsWithOutsidePoints.size() == 0 ||
         facetsWithOutsidePoints.size() == 0 ||
         newFacetsWithOutsidePoints.back()->farthestOutsidePointDistance > facetsWithOutsidePoints.back()->farthestOutsidePointDistance ) {
      facetsWithOutsidePoints.insert( facetsWithOutsidePoints.end(), newFacetsWithOutsidePoints.begin(), newFacetsWithOutsidePoints.end() );
    }
    else {
      // facetsWithOutsidePoints.back() has farther point than newFacetsWithOutsidePoints.back()
      facetsWithOutsidePoints.insert( facetsWithOutsidePoints.begin(), newFacetsWithOutsidePoints.begin(), newFacetsWithOutsidePoints.end() );
    }
  }
  for ( vector< FacetIt >::const_iterator vIt = visibleFacets.begin(); vIt != visibleFacets.end(); ++vIt ) {
    facets.erase( *vIt );
  }

  // Assert that the facets constitute a convex polytope
  updateFacetCenterPoints_( points, facets );
  throwExceptionIfNotConvexPolytope_( facets );
}

/* Anonymous namespace functions */
namespace {
vector< vector< size_t > > computeConvexHull_( const vector< vector< double > >& unperturbedPoints,
                                               double perturbation,
                                               size_t depth )
{
  // Check that the input data is correct
  throwExceptionIfInvalidPerturbation_( perturbation, unperturbedPoints );

  // Perform perturbation of the input points
  const vector< vector< double > >& points = perturbation > 0.0 ? perturbPoints_( unperturbedPoints, perturbation )
                                                                : unperturbedPoints;

  // Get the indices of the extreme points
  const size_t dimension = unperturbedPoints.front().size();
  const vector< size_t > initialPointIndices = getInitialPolytopeVertexIndices_( points );
  assert( initialPointIndices.size() >= dimension + 1 );

  // Create a vector of the extreme points
  vector< vector< double > > initialPoints;
  initialPoints.reserve( initialPointIndices.size() );
  for ( size_t i = 0; i < initialPointIndices.size(); ++i ) {
    initialPoints.push_back( points[ initialPointIndices[ i ] ] );
  }

  // Create simplex using the extreme points
  list< Facet > facets;
  getInitialSimplex_( initialPoints, facets );

  try {
    // Compute the convex hull for the set of all points using the seed polytope
    growConvexHull( initialPoints, facets );
  }
  catch ( invalid_argument e ) {
    // Failed to grow the convex hull. Change perturbation and retry
    cerr << "1 " << e.what() << endl;
    double newPerturbation = perturbation == 0.0 ? 1e-9 : 100 * perturbation;
    throwExceptionIfInvalidPerturbation_( newPerturbation, unperturbedPoints );
    if ( depth > 2 ) {
      throw invalid_argument( "No solution was found although perturbation was increased 3 times." );
    }
    return computeConvexHull_( unperturbedPoints, newPerturbation, depth + 1 );
  }

  // Set the extreme point indices to the indices in the full set of points
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    for ( size_t i = 0; i < fIt->vertexIndices.size(); ++i ) {
      fIt->vertexIndices[ i ] = initialPointIndices[ fIt->vertexIndices[ i ] ];
    }
  }

  try {
    // Compute the convex hull for the set of all points using the seed polytope
    growConvexHull( points, facets );
  }
  catch ( invalid_argument e ) {
    // Failed to grow the convex hull. Change perturbation and retry
    cerr << "2 " << e.what() << endl;
    double newPerturbation = perturbation == 0.0 ? 1e-9 : 100 * perturbation;
    throwExceptionIfInvalidPerturbation_( newPerturbation, unperturbedPoints );
    if ( depth > 2 ) {
      throw invalid_argument( "No solution was found although perturbation was increased 3 times." );
    }
    return computeConvexHull_( unperturbedPoints, newPerturbation, depth + 1 );
  }

  // Assert that the facets' sets of vertices have the correct dimension
  throwExceptionIfNotAllFacetsFullDimensional_( facets, dimension );

  // Construct vector of vertex indices for the facets of the convex hull
  vector< vector< size_t > > vertexIndices( facets.size(), vector< size_t >( dimension ) );
  size_t fi = 0;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt, ++fi ) {
    vertexIndices[ fi ].swap( fIt->vertexIndices );
  }
  return vertexIndices;
}

vector< size_t > getInitialPolytopeVertexIndices_( const vector< vector< double > >& points )
{
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  set< size_t > startPointIndexSet;
  for ( size_t d = 0; d < dimension; ++d ) {
    startPointIndexSet.insert( min_element( points.begin(), points.end(), DimensionComparator( d ) ) - points.begin() );
    startPointIndexSet.insert( max_element( points.begin(), points.end(), DimensionComparator( d ) ) - points.begin() );
  }
  vector< size_t > startPointIndices;
  startPointIndices.reserve( dimension + 1 );
  startPointIndices.insert( startPointIndices.end(), startPointIndexSet.begin(), startPointIndexSet.end() );
  size_t i = 0;
  while ( startPointIndices.size() < dimension + 1 && i < points.size() ) {
    if ( find( startPointIndices.begin(), startPointIndices.end(), i ) == startPointIndices.end() ) {
      startPointIndices.push_back( i );
    }
    ++i;
  }
  return startPointIndices;
}

void getInitialSimplex_( const vector< vector< double > >& points, list< Facet >& facets )
{
  facets.clear();
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  if ( points.size() <= dimension ) {
    throw invalid_argument( "Too few input points to construct convex hull." );
  }

  // Create initial simplex using the (dimension + 1) first points.
  // The facets have vertices [0, ..., dimension - 1], [1, ..., dimension], ..., [dimension, 0, ..., dimension - 2]
  vector< size_t > vertexIndices( dimension );
  for ( size_t i = 0; i <= dimension; ++i ) {
    for ( size_t j = 0; j < dimension; ++j ) {
      vertexIndices[ j ] = ( i + j ) % ( dimension + 1 );
    }
    facets.push_back( Facet( vertexIndices ) );
  }

  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    fIt->isNewFacet = false;
  }

  // Update the facets' neighbors
  for ( FacetIt fIt1 = facets.begin(); fIt1 != facets.end(); ++fIt1 ) {
    for ( FacetIt fIt2 = facets.begin(); fIt2 != facets.end(); ++fIt2 ) {
      if ( fIt1 != fIt2 ) {
        fIt1->neighbors.push_back( fIt2 );
      }
    }
  }
}

vector< vector< double > > perturbPoints_( const vector< vector< double > >& points,
                                           double perturbation )
{
  if ( points.size() == 0 || perturbation == 0.0 ) {
    return points;
  }
  const size_t dimension = points.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( points, dimension );

  vector< vector< double > > perturbedPoints( points.size(), vector< double >( dimension, 0.0 ) );
  for ( size_t pi = 0; pi < points.size(); ++pi ) {
    for ( size_t d = 0; d < dimension; ++d ) {
      perturbedPoints[ pi ][ d ] = points[ pi ][ d ] + ( perturbation * rand() ) / RAND_MAX;
    }
  }
  return perturbedPoints;
}

void updateFacetCenterPoints_( const vector< vector< double > >& points,
                               list< Facet >& facets )
{
  assert( points.size() > 0 );
  const size_t dimension = points.front().size();
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    vector< double > center( dimension, 0.0 );
    for ( vector< size_t >::const_iterator vIt = fIt->vertexIndices.begin(); vIt != fIt->vertexIndices.end(); ++vIt ) {
      const vector< double >& point = points[ *vIt ];
      for ( size_t d = 0; d < dimension; ++d ) {
        center[ d ] += point[ d ];
      }
    }
    for ( size_t d = 0; d < dimension; ++d ) {
      center[ d ] /= fIt->vertexIndices.size();
    }
    fIt->center = center;
  }
}

vector< double > computeOrigin_( const list< Facet >& facets )
{
  assert( facets.size() > 0 );
  const size_t dimension = facets.front().vertexIndices.size();
  vector< double > origin( dimension, 0.0 );
  for ( FacetConstIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    for ( size_t i = 0; i < dimension; ++i ) {
      origin[ i ] += fIt->center[ i ];
    }
  }
  for ( size_t i = 0; i < dimension; ++i ) {
    origin[ i ] /= facets.size();
  }
  return origin;
}

template< template< class T, class All = std::allocator< T > > class Container, class U >
void updateFacetNormalAndOffset_( const vector< vector< double > >& points,
                                  const vector< double >& origin,
                                  Container< U >& facets,
                                  vector< vector< double > >& preallocatedA )
{
  assert( points.size() > 0 );
  const size_t dimension = points.front().size();

  vector< vector< double > >& A = preallocatedA;
  A.resize( dimension );
  for ( size_t i = 0; i + 1 < dimension; ++i ) {
    A[ i ].resize( dimension );
  }

  for ( typename Container< U >::iterator fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    Facet& facet = FacetExtractor< U >::getFacet( fIt );
    assert( facet.vertexIndices.size() == dimension );
    const vector< double >& firstPoint = points[ facet.vertexIndices.front() ];

    for ( size_t i = 0; i + 1 < dimension; ++i ) {
      for ( size_t j = 0; j < dimension; ++j ) {
        A[ i ][ j ] = points[ facet.vertexIndices[ i + 1 ] ][ j ] - firstPoint[ j ];
      }
    }
    vector< double > b( dimension, 0.0 );
    b.back() = 1.0;
    A.back() = b;

    // Solve A x = b
    overwritingSolveLinearSystemOfEquations_( A, b );
    double absSum = 0.0;
    for ( size_t i = 0; i < b.size(); ++i ) {
      absSum += fabs( b[ i ] );
    }

    for ( size_t i = 0; i < b.size(); ++i ) {
      b[ i ] /= absSum;
    }

    facet.normal.swap( b );
    facet.offset = scalarProduct_( facet.normal, points[ facet.vertexIndices.front() ] );

    // Orient normal inwards
    if ( isFacetVisibleFromPoint_( facet, origin ) ) {
      for ( vector< double >::iterator nIt = facet.normal.begin(); nIt != facet.normal.end(); ++nIt ) {
        *nIt = -( *nIt );
      }
      facet.offset = -facet.offset;
    }
  }
}

void initializeOutsideSets_( const vector< vector< double > >& points,
                             list< Facet >& facets )
{
  set< size_t > vertexIndices;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    vertexIndices.insert( fIt->vertexIndices.begin(), fIt->vertexIndices.end() );
  }
  for ( size_t pi = 0; pi < points.size(); ++pi ) {
    if ( vertexIndices.find( pi ) == vertexIndices.end() ) {
      const vector< double >& point = points[ pi ];
      FacetIt farthestFacetIt;
      double maxDistance = 0.0;
      for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
        double distance = distance_( *fIt, point );
        if ( distance > maxDistance ) {
          maxDistance = distance;
          farthestFacetIt = fIt;
        }
      }
      if ( maxDistance > 0.0 ) {
        farthestFacetIt->outsideIndices.push_back( pi );
        if ( maxDistance > farthestFacetIt->farthestOutsidePointDistance ) {
          farthestFacetIt->farthestOutsidePointDistance = maxDistance;
          farthestFacetIt->farthestOutsidePointIndex = pi;
        }
      }
    }
  }
}

size_t getAndEraseFarthestPointFromOutsideSet_( const vector< vector< double > >& points,
                                                Facet& facet )
{
  assert( facet.outsideIndices.size() > 0 );
  vector< size_t >::iterator farthestPointIndexIt =
    find( facet.outsideIndices.begin(), facet.outsideIndices.end(), facet.farthestOutsidePointIndex );
  assert( farthestPointIndexIt != facet.outsideIndices.end() );
  facet.outsideIndices.erase( farthestPointIndexIt );
  return facet.farthestOutsidePointIndex;
}

void getVisibleFacets_( const vector< double >& apex,
                        FacetIt facetIt,
                        vector< FacetIt >& visibleFacets,
                        vector< pair< FacetIt, FacetIt > >& horizon )
{
  facetIt->visible = true;
  facetIt->visitIndex = 0;
  size_t startIndex = visibleFacets.size();
  visibleFacets.push_back( facetIt );
  for ( size_t vi = startIndex; vi < visibleFacets.size(); ++vi ) {
    FacetIt visibleFacetIt = visibleFacets[ vi ];
    const Facet& visibleFacet = *visibleFacetIt;
    for ( size_t ni = 0; ni < visibleFacet.neighbors.size(); ++ni ) {
      FacetIt neighborIt = visibleFacet.neighbors[ ni ];
      Facet& neighbor = *neighborIt;

      if ( neighbor.visitIndex != 0 ) {
        if ( isFacetVisibleFromPoint_( neighbor, apex ) ) {
          visibleFacets.push_back( neighborIt );
          neighbor.visible = true;
        }
      }

      if ( !neighbor.visible ) {
        horizon.push_back( make_pair( visibleFacetIt, neighborIt ) );
      }
      neighbor.visitIndex = 0;
    }
  }
}

size_t getHashValue_( const vector< size_t >& v )
{

  size_t hash = 0;
  for ( size_t i = 0; i < v.size(); ++i ) {
    size_t i2 = ( i + 1 ) * ( i + 1 );
    size_t i4 = i2 * i2;
    size_t i8 = i4 * i4;
    hash += v[ v.size() - i - 1 ] * i8 * i4;
  }
  return hash;
}

void createNewFacets_( size_t apexIndex,
                       const vector< pair< FacetIt, FacetIt > >& horizon,
                       list< Facet >& facets,
                       vector< FacetIt >& visibleFacets,
                       vector< FacetIt >& newFacets,
                       vector< vector< size_t > >& preallocatedPeaks )
{
  assert( horizon.size() > 0 );

  // Construct new facets
  vector< Facet > tmpNewFacets;
  tmpNewFacets.reserve( horizon.size() );

  const size_t dimension = horizon.front().first->vertexIndices.size();
  for ( size_t hi = 0; hi < horizon.size(); ++hi ) {
    const FacetIt visibleFacetIt = horizon[ hi ].first;
    const FacetIt obscuredFacetIt = horizon[ hi ].second;
    assert( visibleFacetIt->visible );
    assert( !obscuredFacetIt->visible );

    // The new facet has the joint vertices of its parent, plus the index of the apex
    assert( find( visibleFacetIt->vertexIndices.begin(), visibleFacetIt->vertexIndices.end(), apexIndex ) == visibleFacetIt->vertexIndices.end() );
    assert( find( obscuredFacetIt->vertexIndices.begin(), obscuredFacetIt->vertexIndices.end(), apexIndex ) == obscuredFacetIt->vertexIndices.end() );
    vector< size_t > vertexIndices( dimension );
    vector< size_t >::const_iterator it =
      set_intersection( visibleFacetIt->vertexIndices.begin(), visibleFacetIt->vertexIndices.end(),
                        obscuredFacetIt->vertexIndices.begin(), obscuredFacetIt->vertexIndices.end(),
                        vertexIndices.begin() );
    assert( size_t( it - vertexIndices.begin() + 1 ) == dimension );
    vertexIndices.back() = apexIndex;
    sort( vertexIndices.begin(), vertexIndices.end() );
    tmpNewFacets.push_back( Facet( vertexIndices ) );
  }

  newFacets.clear();
  newFacets.reserve( horizon.size() );
  for ( size_t hi = 0; hi < horizon.size(); ++hi ) {
    FacetIt newFacetIt;
    if ( visibleFacets.size() > 0 ) {
      // Reuse the space of visible facets, which are to be removed
      *visibleFacets.back() = tmpNewFacets[ hi ];
      newFacetIt = visibleFacets.back();
      visibleFacets.pop_back();
    }
    else {
      facets.push_back( tmpNewFacets[ hi ] );
      newFacetIt = facets.end();
      --newFacetIt;
    }
    // The new facet is neighbor to its obscured parent, and vice versa
    const FacetIt visibleFacetIt = horizon[ hi ].first;
    FacetIt obscuredFacetIt = horizon[ hi ].second;
    vector< FacetIt >::iterator fItIt = find( obscuredFacetIt->neighbors.begin(), obscuredFacetIt->neighbors.end(), visibleFacetIt );
    assert( fItIt != obscuredFacetIt->neighbors.end() );
    *fItIt = newFacetIt;
    obscuredFacetIt->visitIndex = (size_t)-1;
    newFacetIt->neighbors.push_back( obscuredFacetIt );
    newFacets.push_back( newFacetIt );
  }

  vector< vector< size_t > >& peaks = preallocatedPeaks;
  const size_t numOfPeaks = horizon.size() * ( dimension - 1 );
  peaks.resize( numOfPeaks, vector< size_t >( dimension > 1 ? dimension - 2 : 0 ) );
  vector< pair< size_t, pair< FacetIt, vector< size_t >* > > > peakHashes;
  peakHashes.reserve( numOfPeaks );

  size_t peakIndex = 0;
  for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
    FacetIt& newFacetIt = newFacets[ ni ];
    // Find peaks
    for ( size_t i = 0; i < dimension; ++i ) {
      if ( newFacetIt->vertexIndices[ i ] == apexIndex ) {
        continue; // Skip ridges not containing the apex
      }
      size_t pi = 0;
      for ( size_t j = 0; j < dimension; ++j ) {
        // Skip apex
        if ( j == i || newFacetIt->vertexIndices[ j ] == apexIndex ) {
          continue;
        }
        peaks[ peakIndex ][ pi ] = newFacetIt->vertexIndices[ j ];
        ++pi;
      }
      // The vertexIndices are already sorted, so no need to sort them here.
      // If the algorithm is changed to use non-sorted vertices, add the following line:
      // sort( peaks[ peakIndex ].begin(), peaks[ peakIndex ].end() );
      peakHashes.push_back( make_pair( getHashValue_( peaks[ peakIndex ] ), make_pair( newFacetIt, &peaks[ peakIndex ] ) ) );
      ++peakIndex;
    }
  }
  sort( peakHashes.begin(), peakHashes.end(), FirstComparator() );

  // If more than two peaks have the same hash values, it is necessary to sort them by the full vectors
  for ( size_t i = 0; i < peakHashes.size(); ) {
    size_t j = i + 1;
    while ( j < peakHashes.size() && peakHashes[ i ].first == peakHashes[ j ].first ) {
      ++j;
    }
    if ( j > i + 2 ) {
      // More than two peaks have the same hash value
      sort( peakHashes.begin() + i, peakHashes.begin() + j, SecondSecondPtrComparator() );
    }
    i = j;
  }

  // Update neighbors
  for ( size_t ri = 0; ri + 1 < peakHashes.size(); ri += 2 ) {
    FacetIt firstFacetIt = peakHashes[ ri ].second.first;
    FacetIt secondFacetIt = peakHashes[ ri + 1 ].second.first;
    size_t hash1 = peakHashes[ ri ].first;
    size_t hash2 = peakHashes[ ri + 1 ].first;
    const vector< size_t >& peak1 = *peakHashes[ ri ].second.second;
    const vector< size_t >& peak2 = *peakHashes[ ri + 1 ].second.second;
    if ( hash1 != hash2 && peak1 != peak2 ) {
      throw invalid_argument( "Peaks must occur in pairs." );
    }
    firstFacetIt->neighbors.push_back( secondFacetIt );
    secondFacetIt->neighbors.push_back( firstFacetIt );
  }

  for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
    assert( newFacets[ ni ]->neighbors.size() == dimension );
  }
}

void updateOutsideSets_( const vector< vector< double > >& points,
                         const vector< size_t >& outsideIndices,
                         vector< FacetIt >& newFacets )
{
  assert( newFacets.size() > 0 );
  for ( size_t pi = 0; pi < outsideIndices.size(); ++pi ) {
    const size_t pointIndex = outsideIndices[ pi ];
    const vector< double >& point = points[ pointIndex ];
    for ( size_t fi = 0; fi < newFacets.size(); ++fi ) {
      FacetIt newFacetIt = newFacets[ fi ];
      double bestDistance = distance_( *newFacetIt, point );
      if ( bestDistance > 0.0 ) {
        // Found a facet for which the point is an outside point
        // Recursively check whether its neighbors are even farther
        newFacetIt->visitIndex = fi;
        bool checkNeighbors = true;
        while ( checkNeighbors ) {
          checkNeighbors = false;
          for ( size_t ni = 0; ni < newFacetIt->neighbors.size(); ++ni ) {
            FacetIt neighborIt = newFacetIt->neighbors[ ni ];
            if ( !neighborIt->isNewFacet || neighborIt->visitIndex == fi ) {
              continue;
            }
            neighborIt->visitIndex = fi;
            double distance = distance_( *neighborIt, point );
            if ( distance > bestDistance ) {
              bestDistance = distance;
              newFacetIt = neighborIt;
              checkNeighbors = true;
              break;
            }
          }
        }

        newFacetIt->outsideIndices.push_back( pointIndex );
        if ( bestDistance > newFacetIt->farthestOutsidePointDistance ) {
          newFacetIt->farthestOutsidePointDistance = bestDistance;
          newFacetIt->farthestOutsidePointIndex = pointIndex;
        }
        break;
      }
    }
  }
}

vector< size_t > getOutsidePointIndicesFromFacets_( const vector< FacetIt >& facets, size_t startIndex )
{
  size_t numOfUnassignedPoints = 0;
  for ( size_t i = startIndex; i < facets.size(); ++i ) {
    numOfUnassignedPoints += facets[ i ]->outsideIndices.size();
  }
  vector< size_t > unassignedPointIndices;
  unassignedPointIndices.reserve( numOfUnassignedPoints );
  for ( size_t i = startIndex; i < facets.size(); ++i ) {
    const vector< size_t >& outsideIndices = facets[ i ]->outsideIndices;
    unassignedPointIndices.insert( unassignedPointIndices.end(), outsideIndices.begin(), outsideIndices.end() );
  }
  return unassignedPointIndices;
}

bool isFacetVisibleFromPoint_( const Facet& facet,
                               const vector< double >& point )
{
  // Returns true if the point is contained in the open negative halfspace of the facet
  return scalarProduct_( facet.normal, point ) < facet.offset;
}

double distance_( const Facet& facet,
                  const vector< double >& point )
{
  return facet.offset - scalarProduct_( facet.normal, point );
}

double scalarProduct_( const vector< double >& a,
                       const vector< double >& b )
{
  // The scalar product is called a large number of times.
  // The following assert has therefore been removed for speed:
  // assert( a.size() == b.size() );
  //
  // Computing the dot product with the following unrolling gave a speedup of about 15 %
  // for a 3d case with 1e4 randomly distributed points compared to the default method
  switch ( a.size() ) {
    case 2: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ];
    case 3: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ];
    case 4: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] + a[ 3 ] * b[ 3 ];
    case 5: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] + a[ 3 ] * b[ 3 ] + a[ 4 ] * b[ 4 ];
    case 6: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] + a[ 3 ] * b[ 3 ] + a[ 4 ] * b[ 4 ] + a[ 5 ] * b[ 5 ];
    case 7: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] + a[ 3 ] * b[ 3 ] + a[ 4 ] * b[ 4 ] + a[ 5 ] * b[ 5 ] + a[ 6 ] * b[ 6 ];
    case 8: return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] + a[ 3 ] * b[ 3 ] + a[ 4 ] * b[ 4 ] + a[ 5 ] * b[ 5 ] + a[ 6 ] * b[ 6 ] + a[ 7 ] * b[ 7 ];
    default:
    {
      double sum = 0.0;
      for ( size_t i = 0; i < a.size(); ++i ) {
        sum += a[ i ] * b[ i ];
      }
      return sum;
    }
  }
}

void overwritingSolveLinearSystemOfEquations_( vector< vector< double > >& A,
                                               vector< double >& b )
{
  const size_t n = A.size();
  assert( n > 0 );
  for ( size_t i = 0; i < n; ++i ) {
    assert( A[ i ].size() == n );
  }
  assert( b.size() == n );

  // Outer product LU with partial pivoting
  // See Algorithm 3.4.1 in Golub and Van Loan - Matrix Computations, 4th Edition
  for ( size_t k = 0; k < n; ++k ) {
    // Determine mu with k <= mu < n so abs( A( mu, k ) ) = max( A( k:n-1, k ) )
    size_t mu = k;
    double maxValue = fabs( A[ mu ][ k ] );
    for ( size_t i = k + 1; i < n; ++i ) {
      double value = fabs( A[ i ][ k ] );
      if ( value > maxValue ) {
        maxValue = value;
        mu = i;
      }
    }
    if ( maxValue == 0.0 ) {
      // Matrix is singular
      throw invalid_argument( "Singular matrix 1" );
    }
    if ( k != mu ) {
      A[ k ].swap( A[ mu ] );
    }
    // Here, it is utilized that L is not needed
    // (if L is needed, first divide A[ i ][ k ] by A[ k ][ k ], then subtract A[ i ][ k ] * A[ k ][ j ] from A[ i ][ j ])
    for ( size_t i = k + 1; i < n; ++i ) {
      double factor = A[ i ][ k ] / A[ k ][ k ];
      for ( size_t j = k + 1; j < n; ++j ) {
        A[ i ][ j ] -= factor * A[ k ][ j ];
      }
    }
  }
  // LU factorization completed
  // No need to solve Ly = Pb, because b = [0,...,0,1]^T, so y == Pb

  // Solve Ux = y by row-oriented back substitution
  // See Algorithm 3.1.2 in Golub and Van Loan
  for ( size_t j = n; j > 0; --j ) {
    size_t i = j - 1;
    double sum = inner_product( A[ i ].begin() + j, A[ i ].end(), b.begin() + j, 0.0 );
    if ( A[ i ][ i ] != 0.0 ) {
      b[ i ] = ( b[ i ] - sum ) / A[ i ][ i ];
    }
    else {
      // Matrix is singular
      if ( b[ i ] == sum ) {
        // U(i,i) * x(i) == 0.0 and U(i,i) == 0.0 => x(i) == 0.0 is a solution
        b[ i ] = 0.0;
      }
      else {
        // U(i,i) * x(i) != 0.0 but U(i,i) == 0.0 => no solution
        throw invalid_argument( "Singular matrix 2" );
      }
    }
  }
  // b now contains the solution x to Ax = b
}

void throwExceptionIfNotConvexPolytope_( const list< Facet >& facets )
{
  for ( FacetConstIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    for ( size_t i = 0; i < fIt->neighbors.size(); ++i ) {
      if ( isFacetVisibleFromPoint_( *fIt, fIt->neighbors[ i ]->center ) ) {
        throw invalid_argument( "Not a convex polytope" );
      }
    }
  }
}

void throwExceptionIfNotAllFacetsFullDimensional_( const list< Facet >& facets,
                                                   size_t dimension )
{
  for ( FacetConstIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    if ( fIt->vertexIndices.size() != dimension ) {
      throw invalid_argument( "All facets must be full dimensional." );
    }
  }
}

void throwExceptionIfFacetsUseNonExistingVertices_( const list< Facet >& facets,
                                                    const vector< vector< double > >& points )
{
  for ( FacetConstIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    for ( size_t i = 0; i < fIt->vertexIndices.size(); ++i ) {
      if ( points.size() <= fIt->vertexIndices[ i ] ) {
        throw invalid_argument( "All facets must consist of existing vertices." );
      }
    }
  }
}

void throwExceptionIfNotAllPointsHaveCorrectDimension_( const vector< vector< double > >& points,
                                                        size_t dimension )
{
  for ( size_t pi = 0; pi < points.size(); ++pi ) {
    if ( points[ pi ].size() != dimension ) {
      throw invalid_argument( "All points must have the correct dimension." );
    }
  }
}

void throwExceptionIfTooFewPoints_( const vector< vector< double > >& points )
{
  // At least (dimension + 1) points are needed to create a simplex
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  if ( points.size() <= dimension ) {
    throw invalid_argument( "Too few input points to construct convex hull." );
  }
}

void throwExceptionIfInvalidPerturbation_( double perturbation,
                                           const vector< vector< double > >& points )
{
  if ( perturbation < 0.0 ) {
    throw invalid_argument( "Perturbation must be nonnegative." );
  }
  if ( perturbation > 1e-3 ) {
    throw invalid_argument( "Perturbation too large." );
  }
  /*
  if ( perturbation == 0.0 ) {
    return;
  }

  double maxSquareDistance = 0.0;
  for ( size_t i = 0; i < points.size(); ++i ) {
    for ( size_t j = i + 1; j < points.size(); ++j ) {
      double squareDistance = 0.0;
      for ( size_t d = 0; d < points[ i ].size(); ++d ) {
        double difference = ( points[ i ][ d ] - points[ j ][ d ] );
        squareDistance += difference * difference;
      }
      maxSquareDistance = max( maxSquareDistance, squareDistance );
    }
  }
  if ( perturbation * perturbation > maxSquareDistance ) {
    throw invalid_argument( "Perturbation is larger than 1 percent of the maximum distance between points." );
  }
  */
}

} // namespace
} // namespace ConvexHull
