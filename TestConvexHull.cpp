#include <algorithm>
#include <assert.h>
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
      testString_( test )
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

        cout << spaces << output << endl;
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

void testSpeedRandom_()
{
  vector< double > data;
  // these 1e6 points take 0.6 s in qhull
  ifstream fileStream( "points3d-1e4.txt" );
  copy( istream_iterator< double >( fileStream ),
        istream_iterator< double >(),
        back_inserter< vector< double > >( data ) );

  vector< vector< double > > points;
  points.reserve( 1e4 );
  vector< double > point;
  for ( size_t i = 0; i < data.size(); ++i ) {
    point.push_back( data[ i ] );
    if ( point.size() == 3 ) {
      points.push_back( point );
      point.clear();
    }
  }
  for ( size_t i = 0; i < 1; ++i ) {
    vector< vector< size_t > > facets = computeConvexHull( points );
  }
}

void testSpeedUniform_()
{
  vector< vector< double > > points;
  size_t side = 100;
  points.reserve( side * side );
  for ( size_t i = 0; i < side; ++i ) {
    for ( size_t j = 0; j < side; ++j ) {
      vector< double > point( 3 );
      point[ 0 ] = double( i ) - side / 2;
      point[ 1 ] = double( j ) - side / 2;
      point[ 2 ] = point[ 0 ] * point[ 0 ] + point[ 1 ] * point[ 1 ];
      points.push_back( point );
    }
  }
  vector< vector< size_t > > facets = computeConvexHull( points, 1e-8 );
}

int main( int argc, const char* argv[] )
{
  bool standardTest = argc == 1;

  if ( standardTest ) {
    Test t( "Convex hull suite" );
    t.verify( testConvexHull2D_(), "Test 2D" );
    t.verify( testConvexHull3D_(), "Test 3D" );
    t.endTest( true );
  }
  else {
    if ( atoi( argv[ 1 ] ) == 0 ) {
      testSpeedRandom_();
    }
    else {
      testSpeedUniform_();
    }
  }

}
