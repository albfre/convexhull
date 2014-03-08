#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <iterator>
#include <fstream>
#include <vector>

#include "ConvexHull.h"

using namespace std;
using namespace ConvexHull;

class Test {
  public:
    Test( const string& test ) :
      pass_( true ),
      testEnded_( false ),
      testString_( test ),
      start_( clock() )
    {
      string equal;
      for ( size_t i = 0; i < textLength_; ++i ) {
        equal += "=";
      }
      string spaces;
      while ( 2 * spaces.length() + test.length() < textLength_ ) {
        spaces += " ";
      }

      cout << endl;
      cout << equal << endl;
      cout << spaces << test << endl;
      cout << equal << endl << endl;
    };

    ~Test()
    {
      assert( testEnded_ );
    }

    void verify( bool b, const string& s )
    {
      pass_ &= b;
      const string verify = "Verify: ";
      const string spaces = "        ";
      string left = verify + s;
      string output;
      while ( left.length() > textLength_ ) {
        output += left.substr( 0, textLength_ ) + "\n";
        left = spaces + left.substr( textLength_ );
      }
      string dots;
      for ( size_t i = left.length(); i < textLength_; ++i ) {
        dots += ".";
      }
      string passOrFail = b ? "  [OK]" : "[FAIL]";
      output += left + dots + " " + passOrFail;
      cout << output << endl;
    }

    template < typename T >
    void verify( T a, T b, const string& s )
    {
      stringstream ss;
      ss << s << ". Expected value: " << a << ". Actual value: " << b;
      verify( a == b, ss.str() );
    }

    template < typename T >
    void verify( const vector< T >& a, const vector< T >& b, const string& s )
    {
      stringstream ss;
      ss << s << ". Expected value: { ";
      for ( size_t i = 0; i < a.size(); ++i ) {
        ss << a[ i ] << ( i + 1 < a.size() ? ", " : " " );
      }
      ss << "}. Actual value: { ";
      for ( size_t i = 0; i < b.size(); ++i ) {
        ss << b[ i ] << ( i + 1 < b.size() ? ", " : " " );
      }
      ss << "}.";
      verify( a == b, ss.str() );
    }


    bool endTest( bool printPassOrFail = false )
    {
      string passOrFail = pass_ ? "PASSED" : "FAILED";
      string equal;
      for ( size_t i = 0; i < textLength_; ++i ) {
        equal += "=";
      }
      cout << endl;
      cout << equal << endl;

      if ( printPassOrFail ) {
        string output = "Test " + passOrFail;
        string spaces;

        while ( 2 * spaces.length() + output.length() < textLength_ ) {
          spaces += " ";
        }

        stringstream ss;
        ss << "CPU seconds to run test: " << ( clock() - start_ ) / CLOCKS_PER_SEC;
        string time = ss.str();
        string spacesTime;

        while ( 2 * spacesTime.length() + time.length() < textLength_ ) {
          spacesTime += " ";
        }

        cout << spaces << output << endl;
        cout << spacesTime << time << endl;
        cout << equal << endl << endl << endl;
      }
      cout << endl;
      testEnded_ = true;
      return pass_;
    }

  private:
    bool pass_;
    bool testEnded_;
    const string& testString_;
    const double start_;
    static const size_t textLength_ = 70;
};

void sortMatrix( vector< vector< size_t > >& v )
{
  for ( size_t i = 0; i < v.size(); ++i ) {
    sort( v[ i ].begin(), v[ i ].end() );
  }
  sort( v.begin(), v.end() );
}

bool testConvexHull2D_()
{
  Test t( "Computing 2D convex hull." );

  vector< vector< double > > points( 5, vector< double >( 2 ) );
  points[ 1 ][ 0 ] = 1.0;
  points[ 2 ][ 1 ] = 1.0;
  points[ 3 ][ 0 ] = 0.9;
  points[ 3 ][ 1 ] = 0.9;
  points[ 4 ][ 0 ] = -2.0;
  points[ 4 ][ 1 ] = -2.0;

  size_t expected[ 4 ][ 2 ] = {
                                { 1, 3 },
                                { 1, 4 },
                                { 2, 3 },
                                { 2, 4 }
                              };
  vector< vector< size_t > > expectedResults( 4, vector< size_t >( 2 ) );
  for ( size_t i = 0; i < expectedResults.size(); ++i ) {
    expectedResults[ i ].assign( expected[ i ], expected[ i ] + expectedResults[ i ].size() );
  }
  vector< vector< size_t > > facets = computeConvexHull( points );
  sortMatrix( expectedResults );
  sortMatrix( facets );
  t.verify( expectedResults.size(), facets.size(), "Vectors have equal size" );
  for ( size_t i = 0; i < facets.size(); ++i ) {
    t.verify( expectedResults[ i ], facets[ i ], "Facets have equal vertex indices" );
  }

  return t.endTest();
}

bool testConvexHull3D_()
{
  Test t( "Computing 3D convex hull." );

  vector< vector< double > > points( 5, vector< double >( 3 ) );
  // points 0-3 constitute initial simplex
  points[ 1 ][ 0 ] = 1.0;
  points[ 2 ][ 1 ] = 1.0;
  points[ 3 ][ 2 ] = 1.0;

  points[ 4 ][ 0 ] = 0.5;
  points[ 4 ][ 1 ] = -0.4;
  points[ 4 ][ 2 ] = 0.5;

  size_t expected[ 6 ][ 3 ] = {
                                { 0, 1, 2 },
                                { 0, 2, 3 },
                                { 1, 2, 3 },
                                { 0, 1, 4 },
                                { 0, 3, 4 },
                                { 1, 3, 4 }
                              };
  vector< vector< size_t > > expectedResults( 6, vector< size_t >( 3 ) );
  for ( size_t i = 0; i < expectedResults.size(); ++i ) {
    expectedResults[ i ].assign( expected[ i ], expected[ i ] + expectedResults[ i ].size() );
  }
  vector< vector< size_t > > facets = computeConvexHull( points );
  sortMatrix( expectedResults );
  sortMatrix( facets );
  t.verify( expectedResults.size(), facets.size(), "Vectors have equal size" );
  for ( size_t i = 0; i < facets.size(); ++i ) {
    t.verify( expectedResults[ i ], facets[ i ], "Facets have equal vertex indices" );
  }

  return t.endTest();
}

size_t testSpeedRandom_( size_t numOfPoints, size_t dimension, bool print = false )
{
  srand( 123 );
  vector< vector< double > > points;
  points.reserve( numOfPoints );
  stringstream ss;
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    vector< double > point( dimension );
    for ( size_t j = 0; j < dimension; ++j ) {
      point[ j ] = double( rand() ) / RAND_MAX;
      ss << point[ j ] << " ";
    }
    ss << endl;
    points.push_back( point );
  }
  if ( print ) {
    ofstream file( "points.txt" );
    file << dimension << endl << numOfPoints << endl;
    file << ss.str();
    file.close();
  }
  cout << "Convex hull of " << points.size() << " points in " << dimension << "D." << endl << endl;
  double start = clock();
  vector< vector< size_t > > facets = computeConvexHull( points );
  double stop = clock();
  cout << "Number of facets: " << facets.size() << endl;
  cout << "CPU seconds to compute hull (after input): " << ( stop - start ) / CLOCKS_PER_SEC << endl << endl;
  return facets.size();
}

size_t testSpeedUniform_()
{
  vector< vector< double > > points;
  size_t side = 150;
  size_t dimension = 3;
  points.reserve( side * side );
  stringstream ss;
  for ( size_t i = 0; i < side; ++i ) {
    for ( size_t j = 0; j < side; ++j ) {
      vector< double > point( 3 );
      point[ 0 ] = double( i ) - side / 2;
      point[ 1 ] = double( j ) - side / 2;
      point[ 2 ] = point[ 0 ] * point[ 0 ] + point[ 1 ] * point[ 1 ];
      points.push_back( point );
      ss << point[ 0 ] << " " << point[ 1 ] << " " << point[ 2 ] << endl;
    }
  }
  ofstream file( "points.txt" );
  file << dimension << endl << side * side << endl;
  file << ss.str();
  file.close();
  cout << "Convex hull of " << points.size() << " points in " << dimension << "D." << endl << endl;
  double start = clock();
  vector< vector< size_t > > facets = computeConvexHull( points, 1e-8 );
  double stop = clock();
  cout << "Number of facets: " << facets.size() << endl;
  cout << "CPU seconds to compute hull (after input): " << ( stop - start ) / CLOCKS_PER_SEC << endl << endl;
  return facets.size();
}

bool testConvexHullMultiple_()
{
  Test t( "Computing 1D-11D convex hull." );
  size_t expected[ 11 ] = { 2, 11, 38, 149, 534, 1585, 5596, 15353, 41822, 101718, 276556 };

  bool exceptionCaught = false;
  for ( size_t i = 0; i < 11; ++i ) {
    try {
      t.verify( expected[ i ], testSpeedRandom_( 50, i + 1 ), "Number of facets is correct" );
    }
    catch ( exception e ) {
      exceptionCaught = true;
    }
  }
  t.verify( false, exceptionCaught, "No exception" );

  return t.endTest();
}

int main( int argc, const char* argv[] )
{
  bool standardTest = argc == 1;

  if ( standardTest ) {
    Test t( "Convex hull suite" );
    t.verify( testConvexHull2D_(), "Test 2D" );
    t.verify( testConvexHull3D_(), "Test 3D" );
    t.verify( testConvexHullMultiple_(), "Test 1D-11D" );
    t.endTest( true );
  }
  else {
    if ( atoi( argv[ 1 ] ) == 0 ) {
      size_t numOfPoints = size_t( 1e6 );
      size_t dimension = 3;
      if ( argc > 2 ) {
        int np = atoi( argv[ 2 ] );
        assert( np >= 0 );
        numOfPoints = size_t( np );
      }
      if ( argc > 3 ) {
        int dim = atoi( argv[ 3 ] );
        assert( dim > 0 );
        dimension = size_t( dim );
      }
      bool print = true;
      testSpeedRandom_( numOfPoints, dimension, print );
    }
    else if ( atoi( argv[ 1 ] ) == 1 ) {
      testSpeedUniform_();
    }
    else {
      testSpeedUniform_();
      testSpeedRandom_( size_t( 1e6 ), 3 );
    }
  }
}
