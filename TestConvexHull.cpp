#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <iterator>
#include <fstream>
#include <vector>

#include "ConvexHull.h"

class Test {
  public:
    Test(const std::string& test) :
      pass_(true),
      testEnded_(false),
      testString_(test),
      start_(clock())
    {
      std::string equal;
      for (size_t i = 0; i < textLength_; ++i) {
        equal += "=";
      }
      std::string spaces;
      while (2 * spaces.length() + test.length() < textLength_) {
        spaces += " ";
      }

      std::cout << std::endl;
      std::cout << equal << std::endl;
      std::cout << spaces << test << std::endl;
      std::cout << equal << std::endl << std::endl;
    };

    ~Test()
    {
    }

    void verify(bool b, const std::string& s)
    {
      pass_ &= b;
      const auto verify = "Verify: ";
      const auto spaces = "        ";
      auto left = verify + s;
      std::string output;
      while (left.length() > textLength_) {
        output += left.substr(0, textLength_) + "\n";
        left = spaces + left.substr(textLength_);
      }
      std::string dots;
      for (size_t i = left.length(); i < textLength_; ++i) {
        dots += ".";
      }
      auto passOrFail = b ? "  [OK]" : "[FAIL]";
      output += left + dots + " " + passOrFail;
      std::cout << output << std::endl << std::endl;
    }

    template < typename T >
    void verify(T a, T b, const std::string& s)
    {
      std::stringstream ss;
      ss << s << ". Expected value: " << a << ". Actual value: " << b;
      verify(a == b, ss.str());
    }

    template <typename T>
    void verify(const std::vector<T>& a, const std::vector<T>& b, const std::string& s)
    {
      std::stringstream ss;
      ss << s << ". Expected value: { ";
      for (size_t i = 0; i < a.size(); ++i) {
        ss << a[i] << (i + 1 < a.size() ? ", " : " ");
      }
      ss << "}. Actual value: { ";
      for (size_t i = 0; i < b.size(); ++i) {
        ss << b[i] << (i + 1 < b.size() ? ", " : " ");
      }
      ss << "}.";
      verify(a == b, ss.str());
    }


    bool endTest(bool printPassOrFail = false)
    {
      std::string passOrFail = pass_ ? "PASSED" : "FAILED";
      std::string equal;
      for (size_t i = 0; i < textLength_; ++i) {
        equal += "=";
      }
      std::cout << equal << std::endl;

      if (printPassOrFail) {
        auto output = "Test " + passOrFail;
        std::string spaces;

        while (2 * spaces.length() + output.length() < textLength_) {
          spaces += " ";
        }

        std::stringstream ss;
        ss << "CPU seconds to run test: " << std::setprecision(4) << (clock() - start_) / CLOCKS_PER_SEC;
        auto time = ss.str();
        std::string spacesTime;

        while (2 * spacesTime.length() + time.length() < textLength_) {
          spacesTime += " ";
        }

        std::cout << spaces << output << std::endl;
        std::cout << spacesTime << time << std::endl;
        std::cout << equal << std::endl << std::endl << std::endl;
      }
      std::cout << std::endl;
      testEnded_ = true;
      return pass_;
    }

  private:
    bool pass_;
    bool testEnded_;
    const std::string& testString_;
    const double start_;
    static const size_t textLength_ = 70;
};

std::vector<std::vector<size_t>> computeConvexHull(const std::vector<std::vector<double>>& points, double perturbation = 0.0)
{
  const auto printStats = true;
  const auto hull = ConvexHull(points, perturbation, printStats);
  return hull.getVertexIndices();
}

void sortMatrix(std::vector<std::vector<size_t>>& v)
{
  for (auto& vi : v) {
    std::ranges::sort(vi);
  }
  std::ranges::sort(v);
}

bool testConvexHull2D_()
{
  Test t( "Computing 2D convex hull." );

  std::vector<std::vector<double>> points(5, std::vector<double>(2));
  points[1][0] = 1.0;
  points[2][1] = 1.0;
  points[3][0] = 0.9;
  points[3][1] = 0.9;
  points[4][0] = -2.0;
  points[4][1] = -2.0;

  auto expectedResults = std::vector<std::vector<size_t>>{ { 1, 3 },
                                                           { 1, 4 },
                                                           { 2, 3 },
                                                           { 2, 4 } };
  auto facets = computeConvexHull(points);
  sortMatrix(expectedResults);
  sortMatrix(facets);
  t.verify(expectedResults.size(), facets.size(), "Vectors have equal size");
  for (size_t i = 0; i < facets.size(); ++i) {
    t.verify(expectedResults[i], facets[i], "Facets have equal vertex indices");
  }

  return t.endTest();
}

bool testConvexHull3D_()
{
  Test t( "Computing 3D convex hull." );

  std::vector<std::vector<double>> points(5, std::vector<double>(3));
  // points 0-3 constitute initial simplex
  points[1][0] = 1.0;
  points[2][1] = 1.0;
  points[3][2] = 1.0;

  points[4][0] = 0.5;
  points[4][1] = -0.4;
  points[4][2] = 0.5;

  auto expectedResults = std::vector<std::vector<size_t>>{ { 0, 1, 2 },
                                                           { 0, 2, 3 },
                                                           { 1, 2, 3 },
                                                           { 0, 1, 4 },
                                                           { 0, 3, 4 },
                                                           { 1, 3, 4 } };
  auto facets = computeConvexHull(points);
  sortMatrix(expectedResults);
  sortMatrix(facets);
  t.verify(expectedResults.size(), facets.size(), "Vectors have equal size");
  for (size_t i = 0; i < facets.size(); ++i) {
    t.verify(expectedResults[i], facets[i], "Facets have equal vertex indices");
  }

  return t.endTest();
}

bool testConvexHull4D_()
{
  Test t( "Computing 4D convex hull." );
  std::vector<std::vector<double>> input(10, std::vector<double>(4)); // contains a couple of additional instances of the origin
  input[2][0] = 1.0;
  input[3][1] = 1.0;
  input[4][2] = 1.0;
  input[5][3] = 1.0;
  input[6][0] = 0.5;
  input[6][1] = 0.5;
  input[6][2] = -0.4;
  input[6][3] = 0.5;

  auto expectedFacets = std::vector<std::vector<size_t>>{ { 4, 3, 2, 0 },
                                                          { 5, 4, 2, 0 },
                                                          { 4, 5, 3, 0 },
                                                          { 3, 6, 2, 0 },
                                                          { 6, 5, 2, 0 },
                                                          { 5, 6, 3, 0 },
                                                          { 6, 4, 3, 2 },
                                                          { 6, 5, 4, 2 },
                                                          { 5, 6, 4, 3 } };

  auto facets = computeConvexHull(input);
  sortMatrix(expectedFacets);
  sortMatrix(facets);
  t.verify(expectedFacets.size(), facets.size(), "Vectors have equal size");
  for (size_t i = 0; i < facets.size(); ++i) {
    t.verify(expectedFacets[i], facets[i], "Facets have equal vertex indices");
  }

  return t.endTest();
}

std::vector<std::vector<double>> randomPoints(const size_t numOfPoints, const size_t dimension, const bool print)
{
  srand(123);
  std::vector<std::vector<double>> points;
  points.reserve(numOfPoints);
  std::stringstream ss;
  for (size_t i = 0; i < numOfPoints; ++i) {
    std::vector<double> point(dimension);
    for (size_t j = 0; j < dimension; ++j) {
      point[j] = double(rand()) / RAND_MAX;
      ss << point[j] << " ";
    }
    ss << std::endl;
    points.push_back(point);
  }
  if (print) {
    std::ofstream file("points.txt");
    file << dimension << std::endl << numOfPoints << std::endl;
    file << ss.str();
    file.close();
  }
  return points;
}


std::vector<std::vector<double>> uniformPoints(const size_t side)
{
  std::vector<std::vector<double>> points;
  points.reserve(side * side);
  for (size_t i = 0; i < side; ++i) {
    for (size_t j = 0; j < side; ++j) {
      std::vector<double> point(3);
      point[0] = double(i) - side / 2;
      point[1] = double(j) - side / 2;
      point[2] = point[0] * point[0] + point[1] * point[1];
      points.push_back(point);
    }
  }
  return points;
}

size_t testSpeed(const std::vector<std::vector<double>>& points, const size_t loop)
{
  const auto dimension = points.empty() ? 0 : points[0].size();
  std::cout << "Convex hull of " << points.size() << " points in " << dimension << "D." << std::endl << std::endl;
  const double start = clock();
  std::vector<std::vector<size_t>> facets;
  for (size_t i = 0; i < loop; ++i) {
    facets = computeConvexHull(points, 1e-8);
  }
  const double stop = clock();
  std::cout << "Number of facets: " << facets.size() << std::endl;
  std::cout << "CPU seconds to compute hull (after input): " << std::setprecision(4) << (stop - start) / (loop * CLOCKS_PER_SEC) << std::endl << std::endl;
  return facets.size();
}

size_t testSpeedRandom_(size_t numOfPoints, size_t dimension, size_t loop = 1, bool print = false)
{
  const auto points = randomPoints(numOfPoints, dimension, print);
  return testSpeed(points, loop);
}

size_t testSpeedUniform_(const size_t side = 350, const size_t loop = 1)
{
  const auto points = uniformPoints(side);
  return testSpeed(points, loop);
}

bool testConvexHullMultiple_()
{
  Test t( "Computing 1D-10D convex hull." );
  auto expected = std::vector<size_t>{ 2, 11, 38, 149, 534, 1585, 5596, 15353, 41822, 101718, 276556 };

  bool exceptionCaught = false;
  for (size_t i = 0; i < 10; ++i) {
    try {
      t.verify(expected[i], testSpeedRandom_(50, i + 1), "Number of facets is correct");
    }
    catch ( ... ) {
      exceptionCaught = true;
    }
  }
  t.verify(false, exceptionCaught, "No exception");

  return t.endTest();
}

bool testConvexHullHighDim_()
{
  Test t( "Computing 1D-30D convex hull." );

  bool exceptionCaught = false;
  for (size_t i = 0; i < 30; ++i) {
    try {
      t.verify(true, testSpeedRandom_(i + 5, i + 1) > i + 1, "Number of facets is larger than dimension");
    }
    catch ( ... ) {
      exceptionCaught = true;
    }
  }
  t.verify(false, exceptionCaught, "No exception");

  return t.endTest();
}

bool testDisallowedParameters_()
{
  Test t( "Passing disallowed parameters" );
  {
    std::vector<std::vector<double>> points(10, std::vector<double>(3));
    double perturbation = -1e-10;
    std::cout << "Convex hull of " << points.size() << " points in " << points.front().size() << "D with perturbation " << perturbation << "." << std::endl << std::endl;
    bool exceptionCaught = false;
    try {
      auto facets = computeConvexHull(points, perturbation);
    }
    catch ( ... ) {
      exceptionCaught = true;
    }
    t.verify(true, exceptionCaught, "Caught exception");
  }

  {
    std::vector<std::vector<double>> points(10, std::vector<double>(11));
    double perturbation = 1e-10;
    std::cout << "Convex hull of " << points.size() << " points in " << points.front().size() << "D with perturbation " << perturbation << "." << std::endl << std::endl;
    bool exceptionCaught = false;
    try {
      auto facets = computeConvexHull(points, perturbation);
    }
    catch ( ... ) {
      exceptionCaught = true;
    }
    t.verify(true, exceptionCaught, "Caught exception");
  }

  return t.endTest();
}

bool testSimple_()
{
  Test t( "Testing simple" );
  {
    auto points = std::vector<std::vector<double>>{{0, 0, 0},
                                                   {0, 0, 1}, 
                                                   {0, 1, 0}, 
                                                   {0, 1, 1}, 
                                                   {1, 0, 0}, 
                                                   {1, 0, 1}, 
                                                   {1, 1, 0}, 
                                                   {1, 1, 1}};
    auto facets = computeConvexHull(points, 0.);
    std::stringstream ss;
    for (auto& f : facets) {
      for (auto& vi : f) {
        ss << vi << ", ";
      }
      ss << std::endl;
    }
    std::cout << ss.str();
  }
  return t.endTest();
}

size_t getArg(size_t index, const char* argv[])
{
  const auto x = atoi(argv[index]);
  assert(x >= 0);
  return size_t(x);
}

int main(int argc, const char* argv[])
{
  auto standardTest = argc == 1;
  std::vector<size_t> args(argc);
  for (size_t i = 1; i < size_t(argc); ++i) {
    args[i] = getArg(i, argv);
  }

  if (standardTest) {
    Test t( "Convex hull suite" );
    t.verify(testSimple_(), "Test simple 2D");
    t.verify(testConvexHull2D_(), "Test 2D");
    t.verify(testConvexHull3D_(), "Test 3D");
    t.verify(testConvexHull4D_(), "Test 4D");
    t.verify(testConvexHullMultiple_(), "Test 1D-10D");
    t.verify(testConvexHullHighDim_(), "Test 1D-30D");
    t.verify(testDisallowedParameters_(), "Test disallowed parameters");
    t.endTest(true);
  }
  else {
    if (args[1] == 0) {
      size_t numOfPoints = argc > 2 ? args[2] : size_t(1e6);
      size_t dimension = argc > 3 ? args[3] : 3;
      size_t loop = argc > 4 ? args[4] : 1;
      auto print = false;
      testSpeedRandom_(numOfPoints, dimension, loop, print);
    }
    else if (atoi(argv[1]) == 1) {
      size_t loop = argc > 2 ? args[2] : 1;
      testSpeedUniform_(350, loop);
    }
    else {
      testSpeedUniform_();
      testSpeedRandom_(size_t(1e6), 3);
    }
  }
}
