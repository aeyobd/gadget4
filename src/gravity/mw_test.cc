#define BOOST_TEST_MODULE MilkyWayTest

#include "gadgetconfig.h"
#include "grav_milkyway.h"

#include <iostream>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE( test_square ) {
    namespace tt = boost::test_tools;
    BOOST_TEST(square(0.1) == 0.01, tt::tolerance(1e-8));
    BOOST_TEST(square(-2.) == 4., tt::tolerance(1e-8));
    BOOST_TEST(square(0.) == 0., tt::tolerance(1e-8));
}


BOOST_AUTO_TEST_CASE( test_cube ) {
    namespace tt = boost::test_tools;
    BOOST_TEST(cube(0.1) == 0.001, tt::tolerance(1e-8));
    BOOST_TEST(cube(-2.) == -8., tt::tolerance(1e-8));
    BOOST_TEST(cube(0.) == 0., tt::tolerance(1e-8));
}
