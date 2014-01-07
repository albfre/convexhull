#include <algorithm>
#include <iostream>
#include <sstream>

#include "ConvexHull.h"

using namespace std;
using namespace ConvexHull;

namespace {
  bool pass = true;

  void sortMatrix( vector< vector< size_t > >& v )
  {
    for ( size_t i = 0; i < v.size(); ++i ) {
      sort( v[ i ].begin(), v[ i ].end() );
    }
    sort( v.begin(), v.end() );
  }

  void verify( bool b, const string s )
  {
    pass &= b;
    const string verify = "Verify: ";
    const string spaces = "        ";
    string left = verify + s;
    string output = "";
    while ( left.length() > 70 ) {
      output += left.substr( 0, 70 ) + "\n";
      left = spaces + left.substr( 70 );
    }
    string dots = "";
    for ( size_t i = left.length(); i < 70; ++i ) {
      dots += ".";
    }
    string passOrFail = b ? "[PASS]" : "[FAIL]";
    output += left + dots + " " + passOrFail;
    cout << output << endl;
  }

  template < typename T >
  void verify( T a, T b, const string s )
  {
    stringstream ss;
    ss << s << ". Expected value: " << a << ". Actual value: " << b;
    verify( a == b, ss.str() );
  }

  template < typename T >
  void verify( const vector< T >& a, const vector< T >& b, string s )
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

  void endTest()
  {
    string passOrFail = pass ? "PASSED" : "FAILED";
    if ( pass ) {
      cout << endl;
      cout << "=============" << endl;
      cout << " Test " << passOrFail << endl;
      cout << "=============" << endl << endl;
    }
  }
}

int main( int argc, const char* argv[] )
{
  vector< vector< double > > points( 5, vector< double >( 3 ) );
  // points 0-3 constitute initial simplex
  points[ 1 ][ 0 ] = 1.0;
  points[ 2 ][ 1 ] = 1.0;
  points[ 3 ][ 2 ] = 1.0;

  points[ 4 ][ 0 ] = 0.5;
  points[ 4 ][ 1 ] = -0.4;
  points[ 4 ][ 2 ] = 0.5;

  cout << "========================================================" << endl;
  cout << " Test 1. Computing convex hull on nondegenerate points. " << endl;
  cout << "========================================================" << endl << endl;

  for ( size_t i = 0; i < points.size(); ++i ) {
    cout << "Point " << i << ":\t";
    for ( size_t j = 0; j < points[ i ].size(); ++j ) {
      cout << points[ i ][ j ] << "\t";
    }
    cout << endl;
  }
  cout << endl;

  size_t expected[ 4 ][ 3 ] = {
                                { 0, 1, 2 },
                                { 0, 1, 3 },
                                { 0, 2, 3 },
                                { 1, 2, 3 }
                              };
  vector< vector< size_t > > expectedResults( 4, vector< size_t >( 3 ) );
  for ( size_t i = 0; i < expectedResults.size(); ++i ) {
    expectedResults[ i ].assign( expected[ i ], expected[ i ] + expectedResults[ i ].size() );
  }
  vector< vector< size_t > > facets = computeConvexHull( points );
  sortMatrix( expectedResults );
  sortMatrix( facets );
  verify( expectedResults.size(), facets.size(), "Vectors have equal size" );
  for ( size_t i = 0; i < facets.size(); ++i ) {
    verify( expectedResults[ i ], facets[ i ], "Facets have equal vertex indices" );
  }
  endTest();
}
