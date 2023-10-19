#define BOOST_TEST_MODULE MilkyWayTest

#include "gadgetconfig.h"

#include "grav_profiles.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../vectorclass/vectorclass.h"

#include <iostream>
#include <random>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#include <boost/test/included/unit_test.hpp>

namespace tt = boost::test_tools;

typedef std::mt19937 RNG;
RNG rng;
std::normal_distribution<double> rand_normal(0, 1);


double randn() {
    return rand_normal(rng);
}



vector<double> rand_vector(double std=1) {
    double x = randn() * std;
    double y = randn() * std;
    double z = randn() * std;
    return vector<double> {x, y, z};
}


vector<double> normalize(vector<double> vec) {
    return 1./r_sphere(vec) * vec;
}

vector<double> rand_unit() {
    return normalize(rand_vector());
}

/*
 * Sanity tests
 */

BOOST_AUTO_TEST_CASE( test_square ) {
    BOOST_TEST(square(0.1) == 0.01, tt::tolerance(1e-8));
    BOOST_TEST(square(-2.) == 4., tt::tolerance(1e-8));
    BOOST_TEST(square(0.) == 0., tt::tolerance(1e-8));
}


BOOST_AUTO_TEST_CASE( test_cube ) {
    BOOST_TEST(cube(0.1) == 0.001, tt::tolerance(1e-8));
    BOOST_TEST(cube(-2.) == -8., tt::tolerance(1e-8));
    BOOST_TEST(cube(0.) == 0., tt::tolerance(1e-8));
}

BOOST_AUTO_TEST_CASE( test_norm2 ) {
    BOOST_TEST(norm(0., 0.) == 0, tt::tolerance(1e-8));
    BOOST_TEST(norm(3., -4.) == 5., tt::tolerance(1e-8));
    BOOST_TEST(norm(1., 1.) == sqrt(2.), tt::tolerance(1e-8));
    BOOST_TEST(norm(-5., 0.) == 5., tt::tolerance(1e-8));
}



BOOST_AUTO_TEST_CASE( test_origin_acceleration) {
    vector<double> pos = {0, 0, 0};

    double r = r_sphere(NFWAcceleration(pos, 1, 10));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));

    r = r_sphere(DiskAcceleration(pos, 1, 0.5));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));
    r = r_sphere(HernquistAcceleration(pos, 1));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));
}



void test_same_direction(vector<double> a, vector<double>b, double tol=1e-8) {
    vector<double> a_norm = normalize(a);
    vector<double> b_norm = normalize(b);
    double diff = 0;
    for (int i=0; i<3; i++) {
        diff += abs(a_norm[i] - b_norm[i]);
    }
    BOOST_TEST(diff < 1e-8);
}

BOOST_AUTO_TEST_CASE( test_nfw_acceleration_direction ) {
    rng.seed(28);

    vector<double> pos;
    vector<double> acc;

    for (int i=1; i<100; i++) {
        pos = rand_vector(10);
        acc = NFWAcceleration(pos, 1, 10);
        test_same_direction(acc, -1 * pos);
    }
}

BOOST_AUTO_TEST_CASE( test_simple_hernquist ) {
    rng.seed(28);
    double r = 5 + randn();
    vector<double> pos = r * rand_unit();
    double a = r;

    double phi_exp = -1. / (2*a);
    double phi = HernquistPotential(pos, a);
    BOOST_TEST(phi==phi_exp, tt::tolerance(1e-6));
    
    double acc_exp = 1 / square(r + a);
    double acc = r_sphere(HernquistAcceleration(pos, a));
    BOOST_TEST(acc==acc_exp, tt::tolerance(1e-6));
}

// High level mathematical verifications

struct nfw_params {double R; double c;};

double f_nfw(double r, void * params) {
    struct nfw_params * p = (struct nfw_params *)params;
    double R = (p->R);
    double c = (p->c);
    return NFWPotential(r, R, c);
}

BOOST_AUTO_TEST_CASE( test_nfw_acc ) {

    gsl_function F;
    double result, abserr;
    F.function = &f_nfw;
    struct nfw_params p = {2.0, 9.0};
    F.params = &p;

    double r = 1;

    gsl_deriv_forward(&F, r, 1e-8, &result, &abserr);

    double actual = - r * NFWScalarAcceleration(r, p.R, p.c);
    BOOST_TEST(result == actual, tt::tolerance(1e-6));
}



