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
  visitIndex( (size_t)-1 ),
  farthestOutsidePointDistance( 0.0 ),
  visible( false ),
  isNewFacet( true )
{
  neighbors.reserve( vertexIndices.size() );
}

Facet::Facet( vector< size_t >&& _vertexIndices ) :
  vertexIndices( move( _vertexIndices ) ),
  visitIndex( (size_t)-1 ),
  farthestOutsidePointDistance( 0.0 ),
  visible( false ),
  isNewFacet( true )
{
  neighbors.reserve( vertexIndices.size() );
}

typedef list< Facet >::iterator FacetIt;
typedef list< Facet >::const_iterator FacetConstIt;

/* Anonymous namespace function declarations */
namespace {
  int distanceTests;
  int hyperPlanes;
  vector< size_t > INLINE_ATTRIBUTE getInitialPointIndices_( const vector< vector< double > >& points );
  void INLINE_ATTRIBUTE getInitialSimplex_( const vector< vector< double > >& points, const vector< size_t >& startPointIndices, list< Facet >& facets );
  void INLINE_ATTRIBUTE updateFacetCenterPoints_( const vector< vector< double > >& points, list< Facet >& facets );
  vector< double > INLINE_ATTRIBUTE computeOrigin_( const list< Facet >& facets );
  template< template< class T, class All = std::allocator< T > > class Container, class U >
  void INLINE_ATTRIBUTE updateFacetNormalAndOffset_( const vector< vector< double > >& points, const vector< double >& origin, Container< U >& facets, vector< vector< double > >& preallocatedA );
  void INLINE_ATTRIBUTE initializeOutsideSets_( const vector< vector< double > >& points, list< Facet >& facets );
  size_t INLINE_ATTRIBUTE getAndEraseFarthestPointFromOutsideSet_( const vector< vector< double > >& points, Facet& facet );
  void INLINE_ATTRIBUTE getVisibleFacets_( const vector< double >& apex, FacetIt facetIt, vector< FacetIt >& visibleFacets, vector< pair< FacetIt, FacetIt > >& horizon );
  void INLINE_ATTRIBUTE createNewFacets_( size_t apexIndex, const vector< pair< FacetIt, FacetIt > >& horizon, list< Facet >& facets, vector< FacetIt >& visibleFacets, vector< FacetIt >& newFacets, vector< vector< size_t > >& preallocatedPeaks );
  void INLINE_ATTRIBUTE updateOutsideSets_( const vector< vector< double > >& points, const vector< vector< size_t > >& unassignedPointIndices, vector< FacetIt >& newFacets );
  vector< vector< size_t > > INLINE_ATTRIBUTE getOutsidePointIndicesFromFacets_( const vector< FacetIt >& facets, size_t startIndex );
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

  Facet& getFacet_( Facet& f ) { return f; }
  Facet& getFacet_( FacetIt& fIt ) { return *fIt; }

  struct FirstComparator {
    template< class T, class U >
    bool operator() ( const pair< T, U >& p1,
                      const pair< T, U >& p2 ) const { return p1.first < p2.first; } };

  struct FirstSecondSecondPtrComparator {
    template< class T, class U, class V >
    bool operator() ( const pair< T, pair< U, V* > >& p1,
                      const pair< T, pair< U, V* > >& p2 ) const {
                      if ( p1.first != p2.first ) { return p1.first < p2.first; }
                      return *p1.second.second < *p2.second.second; } };

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

  template< size_t N, size_t P > struct Power {
    static const size_t value = N * Power< N, P - 1 >::value;
  };

  template< size_t N > struct Power< N, 0 > {
    static const size_t value = 1;
  };
}

vector< vector< size_t > > computeConvexHull( const vector< vector< double > >& unperturbedPoints,
                                              double perturbation )
{
  distanceTests = 0;
  hyperPlanes = 0;
  // Check that the input data is correct
  throwExceptionIfInvalidPerturbation_( perturbation, unperturbedPoints );
  throwExceptionIfTooFewPoints_( unperturbedPoints );
  const size_t dimension = unperturbedPoints.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( unperturbedPoints, dimension );

  // Set seed to ensure deterministic behavior,
  // but use different seeds for different perturbations
  double* perturbationPtr = &perturbation;
  unsigned int seed = *reinterpret_cast< unsigned int* >( perturbationPtr );
  srand( seed );

  // Perform perturbation of the input points
  struct Point {
    Point() {}
    Point( const vector< double >& coordinates, size_t index ) :
    coordinates_( &coordinates ),
    index_( index )
    {}
    bool operator<( const Point& other ) const {
      return *( this->coordinates_ ) < *( other.coordinates_ );
    }

    bool operator==( const Point& other ) const {
      return *( this->coordinates_ ) == *( other.coordinates_ );
    }

    const vector< double >& getCoordinates() const { return *coordinates_; }
    size_t getIndex() const { return index_; }

    private:
    const vector< double >* coordinates_;
    size_t index_;
  };

  vector< vector< double > > points;
  vector< Point > uniquePoints( unperturbedPoints.size() );
  {
    size_t index = 0;
    transform( unperturbedPoints.cbegin(), unperturbedPoints.cend(), uniquePoints.begin(), [&index] ( const vector< double >& p ) { return Point( p, index++ ); } );
    sort( uniquePoints.begin(), uniquePoints.end() );
    uniquePoints.erase( unique( uniquePoints.begin(), uniquePoints.end() ), uniquePoints.end() );
    points.resize( uniquePoints.size() );
    transform( uniquePoints.cbegin(), uniquePoints.cend(), points.begin(), [] ( const Point& p ) { return p.getCoordinates(); } );
    for ( auto& v : points ) {
      transform( v.cbegin(), v.cend(), v.begin(), [&] ( double d ) { return d + ( perturbation * rand() ) / RAND_MAX; } );
    }
  }

  // Get the indices of the extreme points
  const vector< size_t >& initialPointIndices = getInitialPointIndices_( points );

  // Create simplex using the extreme points
  list< Facet > facets;
  getInitialSimplex_( points, initialPointIndices, facets );

  try {
    // Compute the convex hull for the set of all points using the seed polytope
    growConvexHull( points, facets );
  }
  catch ( invalid_argument e ) {
    // Failed to grow the convex hull. Change perturbation and retry
    cerr << e.what() << endl;
    double newPerturbation = perturbation == 0.0 ? 1e-9 : 10 * perturbation;
    return computeConvexHull( unperturbedPoints, newPerturbation );
  }

  // Assert that the facets' sets of vertices have the correct dimension
  throwExceptionIfNotAllFacetsFullDimensional_( facets, dimension );

  // Construct vector of vertex indices for the facets of the convex hull
  vector< vector< size_t > > vertexIndices( facets.size(), vector< size_t >( dimension ) );
  size_t fi = 0;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt, ++fi ) {
    vertexIndices[ fi ].swap( fIt->vertexIndices );
    transform( vertexIndices[ fi ].cbegin(), vertexIndices[ fi ].cend(), vertexIndices[ fi ].begin(), [&uniquePoints] ( size_t vi ) {
      return uniquePoints[ vi ].getIndex();
    } );
  }
  cerr << "Number of distance tests: " << distanceTests << endl;
  cerr << "Number of hyperplanes created: " << hyperPlanes << endl;
  return vertexIndices;
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
  for ( Facet& facet : facets ) {
    sort( facet.vertexIndices.begin(), facet.vertexIndices.end() );
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
  vector< FacetIt > newFacets;
  while ( !facetsWithOutsidePoints.empty() ) {
    if ( facetsWithOutsidePoints.back()->outsideIndices.empty() ||
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
    const size_t newVisibleFacetsStartIndex = visibleFacets.size();
    getVisibleFacets_( apex, facetIt, visibleFacets, horizon );

    // Get the outside points from the visible facets
    const vector< vector< size_t > >& unassignedPointIndices = getOutsidePointIndicesFromFacets_( visibleFacets, newVisibleFacetsStartIndex );

    // Create new facets from the apex
    newFacets.clear();
    createNewFacets_( apexIndex, horizon, facets, visibleFacets, newFacets, preallocatedPeaks );

    // Update the facet normals and offsets
    updateFacetNormalAndOffset_( points, origin, newFacets, preallocatedA );

    // Assign the points belonging to visible facets to the newly created facets
    updateOutsideSets_( points, unassignedPointIndices, newFacets );

    // Add the new facets with outside points to the vector of all facets with outside points
    for ( FacetIt newFacet : newFacets ) {
      newFacet->isNewFacet = false;
      newFacet->visitIndex = (size_t)-1;
    }
    newFacets.erase( remove_if( newFacets.begin(), newFacets.end(), [] ( const FacetIt& f ) {
        return f->outsideIndices.empty();
      } ), newFacets.end() );

    sort( newFacets.begin(), newFacets.end(), FarthestPointDistanceComparator() );
    if ( newFacets.empty() ||
         facetsWithOutsidePoints.empty() ||
         newFacets.back()->farthestOutsidePointDistance > facetsWithOutsidePoints.back()->farthestOutsidePointDistance ) {
      facetsWithOutsidePoints.insert( facetsWithOutsidePoints.end(), newFacets.begin(), newFacets.end() );
    }
    else {
      // facetsWithOutsidePoints.back() has farther point than newFacetsWithOutsidePoints.back()
      newFacets.insert( newFacets.end(), facetsWithOutsidePoints.begin(), facetsWithOutsidePoints.end() );
      facetsWithOutsidePoints.swap( newFacets );
    }
  }
  for ( const auto& v : visibleFacets ) {
    facets.erase( v );
  }

  // Assert that the facets constitute a convex polytope
  updateFacetCenterPoints_( points, facets );
  throwExceptionIfNotConvexPolytope_( facets );
}

/* Anonymous namespace functions */
namespace {
vector< size_t > getInitialPointIndices_( const vector< vector< double > >& points )
{
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  set< size_t > startPointIndexSet;
  for ( size_t d = 0; d < dimension; ++d ) {
    startPointIndexSet.insert( min_element( points.begin(), points.end(), DimensionComparator( d ) ) - points.begin() );
    startPointIndexSet.insert( max_element( points.begin(), points.end(), DimensionComparator( d ) ) - points.begin() );
  }
  size_t i = 0;
  while ( startPointIndexSet.size() <= dimension && i < points.size()) {
    startPointIndexSet.insert( i );
    ++i;
  }
  return vector< size_t >( startPointIndexSet.begin(), startPointIndexSet.end() );
}

vector< pair< double, pair< size_t, size_t > > > getPairwiseSquaredDistances_( const vector< size_t >& indices,
                                                                               const vector< vector< double > >& points )
{
  vector< pair< double, pair< size_t, size_t > > > distances;
  distances.reserve( indices.size() * ( indices.size() + 1 ) / 2 );
  for ( size_t i = 0; i < indices.size(); ++i ) {
    size_t si = indices[ i ];
    const vector< double >& pointI = points[ si ];
    for ( size_t j = i + 1; j < indices.size(); ++j ) {
      size_t sj = indices[ j ];
      const vector< double >& pointJ = points[ sj ];
      double distance = 0.0;
      for ( size_t k = 0; k < pointI.size(); ++k ) {
        double difference = pointI[ k ] - pointJ[ k ];
        distance += difference * difference;
      }
      distances.push_back( make_pair( distance, make_pair( si, sj ) ) );
    }
  }
  sort( distances.begin(), distances.end() );
  return distances;
}

void getInitialSimplex_( const vector< vector< double > >& points,
                         const vector< size_t >& startPointIndices,
                         list< Facet >& facets )
{
  facets.clear();
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  if ( points.size() <= dimension ) {
    throw invalid_argument( "Too few input points to construct convex hull." );
  }

  vector< pair< double, pair< size_t, size_t > > > distances = getPairwiseSquaredDistances_( startPointIndices, points );

  vector< size_t > sortedIndices;
  sortedIndices.reserve( dimension + 1 );
  while ( sortedIndices.size() <= dimension ) {
    pair< size_t, size_t > indices = distances.back().second;
    distances.pop_back();
    if ( find( sortedIndices.begin(), sortedIndices.end(), indices.first ) == sortedIndices.end() ) {
      sortedIndices.push_back( indices.first );
    }
    if ( find( sortedIndices.begin(), sortedIndices.end(), indices.second ) == sortedIndices.end() ) {
      sortedIndices.push_back( indices.second );
    }
  }

  // Create initial simplex using the (dimension + 1) first points.
  // The facets have vertices [0, ..., dimension - 1], [1, ..., dimension], ..., [dimension, 0, ..., dimension - 2] in sortedIndices.
  vector< size_t > vertexIndices( dimension );
  for ( size_t i = 0; i <= dimension; ++i ) {
    for ( size_t j = 0; j < dimension; ++j ) {
      vertexIndices[ j ] = sortedIndices[ ( i + j ) % ( dimension + 1 ) ];
    }
    facets.push_back( Facet( vertexIndices ) );
    facets.back().isNewFacet = false;
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

void updateFacetCenterPoints_( const vector< vector< double > >& points,
                               list< Facet >& facets )
{
  assert( points.size() > 0 );
  const size_t dimension = points.front().size();
  const double oneOverDimension = 1.0 / dimension;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    vector< double >& center = fIt->center;
    if ( center.size() == dimension ) {
      continue;
    }
    center.assign( dimension, 0.0 );
    for ( vector< size_t >::const_iterator vIt = fIt->vertexIndices.begin(); vIt != fIt->vertexIndices.end(); ++vIt ) {
      const vector< double >& point = points[ *vIt ];
      for ( size_t d = 0; d < dimension; ++d ) {
        center[ d ] += point[ d ];
      }
    }
    for ( size_t d = 0; d < dimension; ++d ) {
      center[ d ] *= oneOverDimension;
    }
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
  for ( size_t i = 0; i < dimension; ++i ) {
    A[ i ].resize( dimension );
  }
  for ( auto& fIt : facets ) {
    Facet& facet = getFacet_( fIt );

    assert( facet.vertexIndices.size() == dimension );
    const vector< double >& firstPoint = points[ facet.vertexIndices.front() ];

    for ( size_t i = 0; i + 1 < dimension; ++i ) {
      for ( size_t j = 0; j < dimension; ++j ) {
        A[ i ][ j ] = points[ facet.vertexIndices[ i + 1 ] ][ j ] - firstPoint[ j ];
      }
    }
    vector< double >& b = facet.normal;
    b.assign( dimension, 0.0 );
    b.back() = 1.0;
    A.back() = b;

    // Solve A x = b
    overwritingSolveLinearSystemOfEquations_( A, b );
    double absSum = 0.0;
    for ( double bi : b ) {
      absSum += fabs( bi );
    }
    absSum = 1.0 / absSum;

    for ( double& bi : b ) {
      bi *= absSum;
    }

    facet.offset = scalarProduct_( facet.normal, points[ facet.vertexIndices.front() ] );
    hyperPlanes++;

    // Orient normal inwards
    if ( isFacetVisibleFromPoint_( facet, origin ) ) {
      transform( facet.normal.cbegin(), facet.normal.cend(), facet.normal.begin(), [] ( double d ) { return -d; } );
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
        const double distance = distance_( *fIt, point );
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

size_t INLINE_ATTRIBUTE getHashValue_( const vector< size_t >& v )
{
  size_t hash = 0;
  vector< size_t >::const_iterator vIt = v.begin();
  switch ( v.size() ) {
    default:
      for ( size_t i = v.size(); i > 15; --i ) {
        size_t i2 = i * i;
        size_t i4 = i2 * i2;
        size_t i8 = i4 * i4;
        hash += *vIt++ * i8 * i4;
      }
    case 15: hash += *vIt++ * Power< 15, 12 >::value;
    case 14: hash += *vIt++ * Power< 14, 12 >::value;
    case 13: hash += *vIt++ * Power< 13, 12 >::value;
    case 12: hash += *vIt++ * Power< 12, 12 >::value;
    case 11: hash += *vIt++ * Power< 11, 12 >::value;
    case 10: hash += *vIt++ * Power< 10, 12 >::value;
    case 9:  hash += *vIt++ * Power< 9, 12 >::value;
    case 8:  hash += *vIt++ * Power< 8, 12 >::value;
    case 7:  hash += *vIt++ * Power< 7, 12 >::value;
    case 6:  hash += *vIt++ * Power< 6, 12 >::value;
    case 5:  hash += *vIt++ * Power< 5, 12 >::value;
    case 4:  hash += *vIt++ * Power< 4, 12 >::value;
    case 3:  hash += *vIt++ * Power< 3, 12 >::value;
    case 2:  hash += *vIt++ * Power< 2, 12 >::value;
    case 1:  hash += *vIt++;
    case 0: ;
  }
  return hash;
}

void INLINE_ATTRIBUTE prepareNewFacets_( size_t apexIndex,
                                         const vector< pair< FacetIt, FacetIt > >& horizon,
                                         list< Facet >& facets,
                                         vector< FacetIt >& visibleFacets,
                                         vector< FacetIt >& newFacets )
{
  vector< Facet > tmpNewFacets;
  tmpNewFacets.reserve( min( horizon.size(), visibleFacets.size() ) );
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
    if ( hi < visibleFacets.size() ) {
      tmpNewFacets.emplace_back( move( vertexIndices ) );
      newFacets.push_back( *( visibleFacets.end() - hi - 1 ) );
    }
    else {
      facets.emplace_back( move( vertexIndices ) );
      newFacets.push_back( facets.end() );
      --newFacets.back();
    }
  }

  // Reuse the space of visible facets, which are to be removed
  for ( size_t hi = 0; hi < tmpNewFacets.size(); ++hi ) {
    swap( *visibleFacets.back(), tmpNewFacets[ hi ] );
    visibleFacets.pop_back();
  }

  // Connect new facets to their neighbors
  for ( size_t hi = 0; hi < horizon.size(); ++hi ) {
    FacetIt newFacetIt = newFacets[ hi ];
    // The new facet is neighbor to its obscured parent, and vice versa
    const FacetIt visibleFacetIt = horizon[ hi ].first;
    FacetIt obscuredFacetIt = horizon[ hi ].second;
    vector< FacetIt >::iterator fItIt = find( obscuredFacetIt->neighbors.begin(), obscuredFacetIt->neighbors.end(), visibleFacetIt );
    assert( fItIt != obscuredFacetIt->neighbors.end() );
    *fItIt = newFacetIt;
    obscuredFacetIt->visitIndex = (size_t)-1;
    newFacetIt->neighbors.push_back( obscuredFacetIt );
  }
}

void INLINE_ATTRIBUTE connectNeighbors_( size_t apexIndex,
                                         const vector< pair< FacetIt, FacetIt > >& horizon,
                                         list< Facet >& facets,
                                         vector< FacetIt >& visibleFacets,
                                         vector< FacetIt >& newFacets,
                                         vector< vector< size_t > >& preallocatedPeaks,
                                         vector< pair< size_t, pair< FacetIt, vector< size_t >* > > >& peakHashes )
{
  const size_t dimension = horizon.front().first->vertexIndices.size();
  vector< vector< size_t > >& peaks = preallocatedPeaks;
  const size_t numOfPeaks = horizon.size() * ( dimension - 1 );
  peaks.resize( numOfPeaks, vector< size_t >( dimension > 1 ? dimension - 2 : 0 ) );
  peakHashes.clear();
  peakHashes.reserve( numOfPeaks );

  size_t peakIndex = 0;
  for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
    FacetIt& newFacetIt = newFacets[ ni ];
    vector< size_t >::const_iterator itStart = newFacetIt->vertexIndices.begin();
    vector< size_t >::const_iterator itEnd = newFacetIt->vertexIndices.end();
    for ( size_t i : newFacetIt->vertexIndices ) {
      if ( i != apexIndex ) {
        vector< size_t >& peak = peaks[ peakIndex ];
        size_t pi = 0;
        for ( size_t j : newFacetIt->vertexIndices ) {
          if ( i != j && j != apexIndex ) {
            peak[ pi++ ] = j;
          }
        }
        // The vertexIndices are already sorted, so no need to sort them here.
        // If the algorithm is changed to use non-sorted vertices, add the following line:
        // sort( peaks[ peakIndex ].begin(), peaks[ peakIndex ].end() );
        const size_t hashVal = getHashValue_( peak );
        peakHashes.push_back( make_pair( hashVal, make_pair( newFacetIt, &peak ) ) );
        ++peakIndex;
      }
    }
  }
  //sort( peakHashes.begin(), peakHashes.end(), FirstSecondSecondPtrComparator() );
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
    const auto& hash1 = peakHashes[ ri ].first;
    const auto& hash2 = peakHashes[ ri + 1 ].first;
    const vector< size_t >& peak1 = *peakHashes[ ri ].second.second;
    const vector< size_t >& peak2 = *peakHashes[ ri + 1 ].second.second;
    if ( hash1 != hash2 && peak1 != peak2 ) {
      throw invalid_argument( "Peaks must occur in pairs." );
    }
    firstFacetIt->neighbors.push_back( secondFacetIt );
    secondFacetIt->neighbors.push_back( firstFacetIt );
  }
}

void createNewFacets_( size_t apexIndex,
                       const vector< pair< FacetIt, FacetIt > >& horizon,
                       list< Facet >& facets,
                       vector< FacetIt >& visibleFacets,
                       vector< FacetIt >& newFacets,
                       vector< vector< size_t > >& preallocatedPeaks )
{
  assert( !horizon.empty() );
  newFacets.clear();
  newFacets.reserve( horizon.size() );
  const size_t dimension = horizon.front().first->vertexIndices.size();

  // Construct new facets
  prepareNewFacets_( apexIndex, horizon, facets, visibleFacets, newFacets );

  // Defining peakHashes here leads to faster solve times for some problems
  vector< pair< size_t, pair< FacetIt, vector< size_t >* > > > peakHashes;
  connectNeighbors_( apexIndex, horizon, facets, visibleFacets, newFacets, preallocatedPeaks, peakHashes );

  for ( const FacetIt& newFacet : newFacets ) {
    assert( newFacet->neighbors.size() == dimension );
  }
}

Facet* assignPointToFarthestFacet_( Facet* facet,
                                    double bestDistance,
                                    size_t pointIndex,
                                    const vector< double >& point,
                                    size_t visitIndex )
{
  // Found a facet for which the point is an outside point
  // Recursively check whether its neighbors are even farther
  facet->visitIndex = visitIndex;
  bool checkNeighbors = true;
  while ( checkNeighbors ) {
    checkNeighbors = false;
    for ( size_t ni = 0; ni < facet->neighbors.size(); ++ni ) {
      Facet& neighbor = *( facet->neighbors[ ni ] );
      if ( !neighbor.isNewFacet || neighbor.visitIndex == visitIndex ) {
        continue;
      }
      neighbor.visitIndex = visitIndex;
      const double distance = distance_( neighbor, point );
      if ( distance > bestDistance ) {
        bestDistance = distance;
        facet = &neighbor;
        checkNeighbors = true;
        break;
      }
    }
  }

  facet->outsideIndices.push_back( pointIndex );
  if ( bestDistance > facet->farthestOutsidePointDistance ) {
    facet->farthestOutsidePointDistance = bestDistance;
    facet->farthestOutsidePointIndex = pointIndex;
  }
  return facet;
}

void updateOutsideSets_( const vector< vector< double > >& points,
                         const vector< vector< size_t > >& visibleFacetOutsideIndices,
                         vector< FacetIt >& newFacets )
{
  assert( !newFacets.empty() );

  for ( size_t vi = 0; vi < visibleFacetOutsideIndices.size(); ++vi ) {
    Facet* facetOfPreviousPoint = nullptr;
    const vector< size_t >& outsideIndices = visibleFacetOutsideIndices[ vi ];
    for ( size_t pi = 0; pi < outsideIndices.size(); ++pi ) {
      const size_t pointIndex = outsideIndices[ pi ];
      const vector< double >& point = points[ pointIndex ];

      if ( facetOfPreviousPoint != nullptr ) {
        double bestDistance = distance_( *facetOfPreviousPoint, point );
        if ( bestDistance > 0.0 ) {
          facetOfPreviousPoint = assignPointToFarthestFacet_( facetOfPreviousPoint, bestDistance, pointIndex, point, pi );
          continue;
        }
      }

      // If the point was not outside the predicted facets, we have to search through all facets
      for ( size_t fi = 0; fi < newFacets.size(); ++fi ) {
        Facet& newFacet = *newFacets[ fi ];
        if ( ( facetOfPreviousPoint == nullptr && facetOfPreviousPoint == &( *newFacets[ fi ] ) )
             ) {
          continue;
        }
        double bestDistance = distance_( newFacet, point );
        if ( bestDistance > 0.0 ) {
          facetOfPreviousPoint = assignPointToFarthestFacet_( &newFacet, bestDistance, pointIndex, point, pi );
          break;
        }
      }
    }
  }
}

vector< vector< size_t > > getOutsidePointIndicesFromFacets_( const vector< FacetIt >& facets,
                                                              size_t startIndex )
{
  vector< vector< size_t > > unassignedPointIndices( facets.size() - startIndex );
  for ( size_t i = startIndex; i < facets.size(); ++i ) {
    unassignedPointIndices[ i - startIndex ].swap( facets[ i ]->outsideIndices );
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

  distanceTests++;
  assert( a.size() == b.size() );
  return inner_product( a.begin(), a.end(), b.begin(), 0.0 );
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
      throw invalid_argument( "Singular matrix 1" );
    }
    if ( k != mu ) {
      A[ k ].swap( A[ mu ] );
    }
    // Here, it is utilized that L is not needed
    // (if L is needed, first divide A[ i ][ k ] by A[ k ][ k ], then subtract A[ i ][ k ] * A[ k ][ j ] from A[ i ][ j ])
    const double invDiag = 1.0 / A[ k ][ k ];
    //size_t numElements = min( n - k - 1, size_t( 15 ) );
    for ( size_t i = k + 1; i < n; ++i ) {
      const double factor = A[ i ][ k ] * invDiag;
      // The loop unrolling below is equivalent to this:
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
    const size_t i = j - 1;
    const double sum = inner_product( A[ i ].begin() + j, A[ i ].end(), b.begin() + j, 0.0 );

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
  if ( perturbation == 0.0 ) {
    return;
  }

  const vector< size_t >& initialPointIndices = getInitialPointIndices_( points );
  vector< pair< double, pair< size_t, size_t > > > distances = getPairwiseSquaredDistances_( initialPointIndices, points );
  double maxDistance = sqrt( distances.back().first );
  if ( perturbation > 0.01 * maxDistance ) {
    throw invalid_argument( "Perturbation is larger than 1 percent of the maximum distance between points." );
  }
}

} // namespace
} // namespace ConvexHull
