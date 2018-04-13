#include "gtest/gtest.h"

#include "polyATests.hpp"
#include "knownAluTests.hpp"


GTEST_API_ int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
