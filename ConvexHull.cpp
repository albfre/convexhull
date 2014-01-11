/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <set>

/* HEADER */
#include "ConvexHull.h"

using namespace std;

namespace ConvexHull {

Facet::Facet( const vector< size_t >& _vertexIndices ) : vertexIndices( _vertexIndices ), visible( false ), visited( false ) {}
typedef list< Facet >::iterator FacetIt;
typedef list< Facet >::const_iterator FacetConstIt;

/* Anonymous namespace functions declarations */
namespace {
  list< Facet > getInitialSimplex_( const vector< vector< double > >& points );
  vector< vector< double > > perturbPoints_( const vector< vector< double > >& points, double perturbation );
  void updateFacetCenterPoints_( const vector< vector< double > >& points, list< Facet >& facets );
  vector< double > computeOrigin_( const list< Facet >& facets );
  void updateFacetNormalAndOffset_( const vector< vector< double > >& points, const vector< double >& origin, Facet& facet, vector< vector< double > >& preallocatedA );
  void initializeOutsideSets_( const vector< vector< double > >& points, list< Facet >& facets );
  size_t getAndEraseFarthestPointFromOutsideSet_( const vector< vector< double > >& points, Facet& facet );
  vector< FacetIt > getVisibleFacets_( const vector< double >& apex, FacetIt facetIt, vector< pair< FacetIt, FacetIt > >& horizon );
  vector< FacetIt > createNewFacets_( size_t apexIndex, const vector< pair< FacetIt, FacetIt > >& horizon, list< Facet >& facets );
  void updateOutsideSets_( const vector< vector< double > >& points, const vector< FacetIt >& visibleFacets, vector< FacetIt >& newFacets );
  bool isFacetVisibleFromPoint_( const Facet& facet, const vector< double >& point );
  double scalarProduct_( const vector< double >& a, const vector< double >& b );
  void overwritingSolveLinearSystemOfEquations_( vector< vector< double > >& A, vector< double >& b );

  void throwExceptionIfNotConvexPolytope_( const list< Facet >& facets );
  void throwExceptionIfNotAllFacetsFullDimensional_( const list< Facet >& facets, size_t dimension );
  void throwExceptionIfFacetsUseNonExistingVertices_( const list< Facet >& facets, const vector< vector< double > >& points );
  void throwExceptionIfNotAllPointsHaveCorrectDimension_( const vector< vector< double > >& points, size_t dimension );
  void throwExceptionIfTooFewPoints_( const vector< vector< double > >& points );
  void throwExceptionIfInvalidPerturbation_( double perturbation, const vector< vector< double > >& points );
}

vector< vector< size_t > > computeConvexHull( const vector< vector< double > >& unperturbedPoints,
                                              double perturbation )
{
  // Check that the input data is correct
  throwExceptionIfTooFewPoints_( unperturbedPoints );
  const size_t dimension = unperturbedPoints.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( unperturbedPoints, dimension );
  throwExceptionIfInvalidPerturbation_( perturbation, unperturbedPoints );

  // Perform perturbation of the input points
  const vector< vector< double > >& points = perturbation > 0.0 ? perturbPoints_( unperturbedPoints, perturbation )
                                                                : unperturbedPoints;

  // Get the initial simplex to use as a seed polytope
  list< Facet > facets = getInitialSimplex_( points );

  try {
    // Compute the convex hull for the set of all points using the seed polytope
    growConvexHull( points, facets );
  }
  catch ( invalid_argument e ) {
    // Failed to grow the convex hull. Change perturbation and retry
    double newPerturbation = perturbation == 0.0 ? 1e-9 : 100 * perturbation;
    throwExceptionIfInvalidPerturbation_( newPerturbation, unperturbedPoints );
    return computeConvexHull( unperturbedPoints, newPerturbation );
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

void growConvexHull( const vector< vector< double > >& points,
                     list< Facet >& facets )
{
  // Check that the input data is correct
  throwExceptionIfTooFewPoints_( points );
  const size_t dimension = points.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( points, dimension );
  throwExceptionIfNotAllFacetsFullDimensional_( facets, dimension );
  throwExceptionIfFacetsUseNonExistingVertices_( facets, points );

  // Update the facet center points to the mean of the vertex points
  updateFacetCenterPoints_( points, facets );

  // Compute origin as the mean of the center points of the seed facets
  const vector< double > origin = computeOrigin_( facets );

  // Compute inwards-oriented facet normals
  vector< vector< double > > preallocatedA;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    updateFacetNormalAndOffset_( points, origin, *fIt, preallocatedA );
  }
  throwExceptionIfNotConvexPolytope_( facets );

  // Compute outside sets
  initializeOutsideSets_( points, facets );

  // Create a list of all facets that have outside points
  list< FacetIt > facetsWithOutsidePoints;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    if ( fIt->outsideIndices.size() > 0 ) {
      facetsWithOutsidePoints.push_back( fIt );
    }
  }

  // Sort vertex indices
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    sort( fIt->vertexIndices.begin(), fIt->vertexIndices.end() );
  }

  // Create new facets using the outside points
  while ( facetsWithOutsidePoints.size() > 0 ) {
    FacetIt facetIt = facetsWithOutsidePoints.back();
    facetsWithOutsidePoints.pop_back();

    // From the outside set of the current facet, find the farthest point
    const size_t apexIndex = getAndEraseFarthestPointFromOutsideSet_( points, *facetIt );
    const vector< double >& apex = points[ apexIndex ];

    // Find the set of facets that are visible from the point to be added
    vector< pair< FacetIt, FacetIt > > horizon; // visible-invisible neighboring facet pairs
    const vector< FacetIt > visibleFacets = getVisibleFacets_( apex, facetIt, horizon );

    // Create new facets from the apex
    vector< FacetIt > newFacets = createNewFacets_( apexIndex, horizon, facets );

    // Update the facet normals and offsets
    for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
      updateFacetNormalAndOffset_( points, origin, *newFacets[ ni ], preallocatedA );
    }

    // Assign the points belonging to visible facets to the newly create facets
    updateOutsideSets_( points, visibleFacets, newFacets );

    // Remove the visible facets from the sets of facets
    for ( size_t vi = 0; vi < visibleFacets.size(); ++vi ) {
      facetsWithOutsidePoints.remove( visibleFacets[ vi ] );
      facets.erase( visibleFacets[ vi ] );
    }

    // Add the new facets with outside points
    for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
      if ( newFacets[ ni ]->outsideIndices.size() > 0 ) {
        facetsWithOutsidePoints.push_back( newFacets[ ni ] );
      }
    }
  }

  // Assert that the facets constitute a convex polytope
  updateFacetCenterPoints_( points, facets );
  throwExceptionIfNotConvexPolytope_( facets );
}

/* Anonymous namespace functions */
namespace {

list< Facet > getInitialSimplex_( const vector< vector< double > >& points )
{
  const size_t dimension = points.size() == 0 ? 0 : points.front().size();
  if ( points.size() <= dimension ) {
    throw invalid_argument( "Too few input points to construct convex hull." );
  }
  list< Facet > facets;

  // Create initial simplex using the (dimension + 1) first points.
  // The facets have vertices [0, ..., dimension - 1], [1, ..., dimension], ..., [dimension, 0, ..., dimension - 2]
  for ( size_t i = 0; i <= dimension; ++i ) {
    vector< size_t > vertexIndices( dimension );
    for ( size_t j = 0; j < dimension; ++j ) {
      vertexIndices[ j ] = ( i + j ) % ( dimension + 1 );
    }
    facets.push_back( Facet( vertexIndices ) );
  }

  // Update the facets' neighbors
  for ( FacetIt fIt1 = facets.begin(); fIt1 != facets.end(); ++fIt1 ) {
    for ( FacetIt fIt2 = facets.begin(); fIt2 != facets.end(); ++fIt2 ) {
      if ( fIt1 != fIt2 ) {
        fIt1->neighbors.push_back( fIt2 );
      }
    }
  }

  return facets;
}

vector< vector< double > > perturbPoints_( const vector< vector< double > >& points,
                                           double perturbation )
{
  if ( points.size() == 0 || perturbation == 0.0 ) {
    return points;
  }
  const size_t dimension = points.front().size();
  throwExceptionIfNotAllPointsHaveCorrectDimension_( points, dimension );

  // always use the same seed to ensure deterministic behavior
  srand( 4711 );
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

void updateFacetNormalAndOffset_( const vector< vector< double > >& points,
                                  const vector< double >& origin,
                                  Facet& facet,
                                  vector< vector< double > >& preallocatedA )
{
  assert( points.size() > 0 );
  const size_t dimension = points.front().size();
  assert( facet.vertexIndices.size() == dimension );
  vector< vector< double > >& A = preallocatedA;
  A.resize( dimension );
  for ( size_t i = 0; i < dimension; ++i ) {
    A[ i ].resize( dimension );
    for ( size_t j = 0; j < dimension; ++j ) {
      A[ i ][ j ] = points[ facet.vertexIndices[ i ] ][ j ] - origin[ j ];
    }
  }
  vector< double > b( dimension, 1.0 );

  // Solve A x = b
  overwritingSolveLinearSystemOfEquations_( A, b );
  double sum = accumulate( b.begin(), b.end(), 0.0 );
  for ( size_t i = 0; i < b.size(); ++i ) {
    b[ i ] /= sum;
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

void initializeOutsideSets_( const vector< vector< double > >& points,
                             list< Facet >& facets )
{
  vector< size_t > unassignedPointIndices;
  set< size_t > vertexIndices;
  for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
    vertexIndices.insert( fIt->vertexIndices.begin(), fIt->vertexIndices.end() );
  }
  for ( size_t pi = 0; pi < points.size(); ++pi ) {
    if ( vertexIndices.find( pi ) == vertexIndices.end() ) {
      unassignedPointIndices.push_back( pi );
    }
  }

  for ( vector< size_t >::const_iterator pIt = unassignedPointIndices.begin(); pIt != unassignedPointIndices.end(); ++pIt ) {
    for ( FacetIt fIt = facets.begin(); fIt != facets.end(); ++fIt ) {
      if ( isFacetVisibleFromPoint_( *fIt, points[ *pIt ] ) ) {
        fIt->outsideIndices.push_back( *pIt );
        break;
      }
    }
  }
}

size_t getAndEraseFarthestPointFromOutsideSet_( const vector< vector< double > >& points,
                                                Facet& facet )
{
  assert( facet.outsideIndices.size() > 0 );
  vector< size_t >::iterator farthestPointIndexIt = facet.outsideIndices.begin();
  double maxOffset = scalarProduct_( facet.normal, points[ *farthestPointIndexIt ] );
  for ( vector< size_t >::iterator pIt = facet.outsideIndices.begin(); pIt != facet.outsideIndices.end(); ++pIt ) {
    double offset = scalarProduct_( facet.normal, points[ *pIt ] );
    if ( offset > maxOffset ) {
      maxOffset = offset;
      farthestPointIndexIt = pIt;
    }
  }
  size_t index = *farthestPointIndexIt;
  facet.outsideIndices.erase( farthestPointIndexIt );
  return index;
}

vector< FacetIt > getVisibleFacets_( const vector< double >& apex,
                                     FacetIt facetIt,
                                     vector< pair< FacetIt, FacetIt > >& horizon )
{
  facetIt->visible = true;
  facetIt->visited = true;
  vector< FacetIt > visibleFacets( 1, facetIt );
  for ( size_t vi = 0; vi < visibleFacets.size(); ++vi ) {
    FacetIt visibleFacetIt = visibleFacets[ vi ];
    const Facet& visibleFacet = *visibleFacetIt;
    for ( size_t ni = 0; ni < visibleFacet.neighbors.size(); ++ni ) {
      FacetIt neighborIt = visibleFacet.neighbors[ ni ];
      Facet& neighbor = *neighborIt;

      if ( !neighbor.visited ) {
        if ( isFacetVisibleFromPoint_( neighbor, apex ) ) {
          visibleFacets.push_back( neighborIt );
          neighbor.visible = true;
        }
      }

      if ( !neighbor.visible ) {
        horizon.push_back( make_pair( visibleFacetIt, neighborIt ) );
      }
      neighbor.visited = true;
    }
  }
  return visibleFacets;
}

struct SecondComparator { bool operator() ( const pair< FacetIt, vector< size_t > >& p1,
                                            const pair< FacetIt, vector< size_t > >& p2 )
                          { return p1.second < p2.second; } };

vector< FacetIt > createNewFacets_( size_t apexIndex,
                                    const vector< pair< FacetIt, FacetIt > >& horizon,
                                    list< Facet >& facets )
{
  assert( horizon.size() > 0 );

  // Construct new facets
  vector< FacetIt > newFacets;
  vector< pair< FacetIt, vector< size_t > > > facetPeakPairs;
  const size_t dimension = horizon.front().first->vertexIndices.size();
  facetPeakPairs.reserve( horizon.size() * ( dimension + 1 ) );
  for ( size_t hi = 0; hi < horizon.size(); ++hi ) {
    FacetIt visibleFacetIt = horizon[ hi ].first;
    FacetIt obscuredFacetIt = horizon[ hi ].second;
    assert( visibleFacetIt->visible );
    assert( !obscuredFacetIt->visible );

    // The new facet has the joint vertices of its parent, plus the index of the apex
    assert( find( visibleFacetIt->vertexIndices.begin(), visibleFacetIt->vertexIndices.end(), apexIndex ) == visibleFacetIt->vertexIndices.end() );
    assert( find( obscuredFacetIt->vertexIndices.begin(), obscuredFacetIt->vertexIndices.end(), apexIndex ) == obscuredFacetIt->vertexIndices.end() );
    vector< size_t > vertexIndices( dimension );
    vector< size_t >::iterator it =
      set_intersection( visibleFacetIt->vertexIndices.begin(), visibleFacetIt->vertexIndices.end(),
                        obscuredFacetIt->vertexIndices.begin(), obscuredFacetIt->vertexIndices.end(),
                        vertexIndices.begin() );
    assert( size_t( it - vertexIndices.begin() + 1 ) == dimension );
    vertexIndices.back() = apexIndex;
    sort( vertexIndices.begin(), vertexIndices.end() );

    facets.push_back( Facet( vertexIndices ) );
    FacetIt newFacetIt = facets.end();
    --newFacetIt;
    newFacets.push_back( newFacetIt );

    // The new facet is neighbor to its obscured parent, and vice versa
    newFacetIt->neighbors.push_back( obscuredFacetIt );
    vector< FacetIt >::iterator fItIt = find( obscuredFacetIt->neighbors.begin(), obscuredFacetIt->neighbors.end(), visibleFacetIt );
    assert( fItIt != obscuredFacetIt->neighbors.end() );
    *fItIt = newFacetIt;

    obscuredFacetIt->visited = false;

    // Find peaks
    for ( size_t i = 0; i < dimension; ++i ) {
      if ( newFacetIt->vertexIndices[ i ] == apexIndex ) {
        continue; // Skip ridges not containing the apex
      }
      vector< size_t > peak;
      peak.reserve( dimension - 1 );
      for ( size_t j = 0; j < dimension; ++j ) {
        // Skip apex
        if ( j == i || newFacetIt->vertexIndices[ j ] == apexIndex ) {
          continue;
        }
        peak.push_back( newFacetIt->vertexIndices[ j ] );
      }
      sort( peak.begin(), peak.end() );
      facetPeakPairs.push_back( make_pair( newFacetIt, peak ) );
    }
  }
  sort( facetPeakPairs.begin(), facetPeakPairs.end(), SecondComparator() );

  // Update neighbors
  for ( size_t ri = 0; ri + 1 < facetPeakPairs.size(); ri += 2 ) {
    FacetIt firstFacetIt = facetPeakPairs[ ri ].first;
    FacetIt secondFacetIt = facetPeakPairs[ ri + 1 ].first;
    const vector< size_t >& peak1 = facetPeakPairs[ ri ].second;
    const vector< size_t >& peak2 = facetPeakPairs[ ri + 1 ].second;
    if ( peak1 != peak2 ) {
      throw invalid_argument( "Peaks must occur in pairs." );
    }
    firstFacetIt->neighbors.push_back( secondFacetIt );
    secondFacetIt->neighbors.push_back( firstFacetIt );
  }

  for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
    assert( newFacets[ ni ]->neighbors.size() == dimension );
  }

  return newFacets;
}

void updateOutsideSets_( const vector< vector< double > >& points,
                         const vector< FacetIt >& visibleFacets,
                         vector< FacetIt >& newFacets )
{
  size_t numOfUnassignedPoints = 0;
  for ( size_t vi = 0; vi < visibleFacets.size(); ++vi ) {
    numOfUnassignedPoints += visibleFacets[ vi ]->outsideIndices.size();
  }
  vector< size_t > unassignedPointIndices;
  unassignedPointIndices.reserve( numOfUnassignedPoints );
  for ( size_t vi = 0; vi < visibleFacets.size(); ++vi ) {
    const vector< size_t >& outsideIndices = visibleFacets[ vi ]->outsideIndices;
    unassignedPointIndices.insert( unassignedPointIndices.end(), outsideIndices.begin(), outsideIndices.end() );
  }

  for ( vector< size_t >::const_iterator pIt = unassignedPointIndices.begin(); pIt != unassignedPointIndices.end(); ++pIt ) {
    for ( size_t ni = 0; ni < newFacets.size(); ++ni ) {
      if ( isFacetVisibleFromPoint_( *newFacets[ ni ], points[ *pIt ] ) ) {
        newFacets[ ni ]->outsideIndices.push_back( *pIt );
        break;
      }
    }
  }
}

bool isFacetVisibleFromPoint_( const Facet& facet,
                               const vector< double >& point )
{
  // Returns true if the point is contained in the open negative halfspace of the facet
  return scalarProduct_( facet.normal, point ) < facet.offset;
}

double scalarProduct_( const vector< double >& a,
                       const vector< double >& b )
{
  assert( a.size() == b.size() );
  return inner_product( a.begin(), a.end(), b.begin(), 0.0 );
}

void overwritingSolveLinearSystemOfEquations_( vector< vector< double > >& A,
                                               vector< double >& b )
{
  const string singularMatrix = "Singular matrix";
  const size_t n = A.size();
  assert( n > 0 );
  for ( size_t i = 0; i < A.size(); ++i ) {
    assert( A[ i ].size() == n );
  }
  assert( b.size() == n );
  vector< size_t > pivot;
  pivot.resize( n - 1 );

  // Outer product LU with partial pivoting
  // See Algorithm 3.4.1 in Golub and Van Loan - Matrix Computations, 4th Edition
  for ( size_t k = 0; k < n - 1; ++k ) {
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
    pivot[ k ] = mu;
    A[ k ].swap( A[ mu ] );

    if ( A[ k ][ k ] != 0.0 ) {
      // rho = k + 1:n - 1
      // A(rho,k) = A(rho,k) / A(k,k)
      double factor = 1.0 / A[ k ][ k ];
      for ( size_t rho = k + 1; rho < n; ++rho ) {
        A[ rho ][ k ] *= factor;
      }

      // A(rho,rho) = A(rho,rho) - A(rho,k) * A(k,rho)
      for ( size_t i = k + 1; i < n; ++i ) {
        for ( size_t j = k + 1; j < n; ++j ) {
          A[ i ][ j ] -= A[ i ][ k ] * A[ k ][ j ];
        }
      }
    }
    else {
      // Matrix is singular
      throw invalid_argument( singularMatrix );
    }
  }
  if ( A[ n - 1 ][ n - 1 ] == 0.0 ) {
    // Matrix is singular
    throw invalid_argument( singularMatrix );
  }
  // LU factorization completed

  // Replace b by Pb
  for ( size_t j = 0; j < n - 1; ++j ) {
    swap( b[ j ], b[ pivot[ j ] ] );
  }

  // Solve Ly = Pb by column-oriented forward substitution
  // See Algorithm 3.1.3 in Golub and Van Loan (utilizing that L is unit lower triangular)
  for ( size_t j = 0; j < n - 1; ++j ) {
    for ( size_t k = j + 1; k < n; ++k ) {
      b[ k ] -= b[ j ] * A[ k ][ j ];
    }
  }

  // Solve Ux = y by row-oriented back substitution
  // See Algorithm 3.1.2 in Golub and Van Loan
  for ( size_t j = n; j > 0; --j ) {
    size_t i = j - 1;
    double sum = inner_product( A[ i ].begin() + j, A[ i ].end(), b.begin() + j, 0.0 );
    if ( A[ i ][ i ] == 0.0 ) {
      // Matrix is singular
      if ( b[ i ] == sum ) {
        // U(i,i) * x(i) == 0.0 and U(i,i) == 0.0 => x(i) == 0.0 is a solution
        b[ i ] = 0.0;
      }
      else {
        // U(i,i) * x(i) != 0.0 but U(i,i) == 0.0 => no solution
        throw invalid_argument( singularMatrix );
      }
    }
    else {
      b[ i ] = ( b[ i ] - sum ) / A[ i ][ i ];
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
}

} // namespace
} // namespace ConvexHull
