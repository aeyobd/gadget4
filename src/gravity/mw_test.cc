#define BOOST_TEST_MODULE test mw potentials
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "gadgetconfig.h"

#include "grav_profiles.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../vectorclass/vectorclass.h"

#include <iostream>
#include <random>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/unit_test_parameters.hpp>


// declarations
namespace tt = boost::test_tools;
namespace bdata = boost::unit_test::data;

typedef std::mt19937 RNG;
typedef vector<double> fvector;
RNG rng;
std::normal_distribution<double> rand_normal(0, 1);

int Npoints = 100;
double r_min_exp = -3;
double r_max_exp = 4;
double sigma=0.05;

struct LogConfig {
    LogConfig() : test_log("test.log"){
        boost::unit_test::unit_test_log.set_stream(test_log);
        // boost::unit_test::unit_test_log.set_format(HRF);
    }
    ~LogConfig() {
        boost::unit_test::unit_test_log.set_stream(std::cout);
    }
    std::ofstream test_log;
};

BOOST_TEST_GLOBAL_CONFIGURATION(LogConfig);



// HELPER FUNCTIONS


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



// Data makers
//
std::vector<double> test_radii() {
    std::vector<double> radii (Npoints);
    double exponent, r;

    for (int i=0; i<Npoints; i++) {
        exponent = r_min_exp + (r_max_exp-r_min_exp)/Npoints * i;
        r = pow(10., exponent);
        r *= (1 + sigma*randn());
        radii[i] = r;
    };

    return radii;
}


std::vector<double> large_radii() {
    std::vector<double> radii (Npoints);
    double exponent, r;

    for (int i=0; i<Npoints; i++) {
        exponent = r_min_exp + (r_max_exp-r_min_exp)/Npoints * i;
        exponent += 8;
        r = pow(10., exponent);
        r *= (1 + sigma*randn());
        radii[i] = r;
    };

    return radii;
}

auto rand_a = bdata::random(1.0, 10.0);
auto rand_b = bdata::random(0.1, 2.0);
auto rand_c = bdata::random(2.0, 20.0);



std::vector<fvector> vectors() {
    std::vector<fvector> vecs (Npoints);
    std::vector<double> radii = test_radii();

    for(int i=0; i<Npoints; i++) {
        vecs[i] = radii[i] * rand_unit();
    }
    return vecs;
}


// required for data tests
static inline std::ostream& operator<<(std::ostream& os, const fvector& pos) {
    os << "vec";
    return os;
}





/*
 * Inline function tests (just for sanity)
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


BOOST_AUTO_TEST_CASE( test_norm3 ) {
    BOOST_TEST(norm(0., 0., 0.) == 0, tt::tolerance(1e-8));
    BOOST_TEST(norm(3., 0., -4.) == 5., tt::tolerance(1e-8));
    BOOST_TEST(norm(1., -1., 1.) == sqrt(3.), tt::tolerance(1e-8));
    BOOST_TEST(norm(-5., 6., 7.) == sqrt(110.), tt::tolerance(1e-8));
}


BOOST_AUTO_TEST_CASE( test_math ) {
    BOOST_TEST(log(M_E) == 1, tt::tolerance(1e-12));
    BOOST_TEST(pow(10., 3) == 1000., tt::tolerance(1e-12));
}


// ================================
// special point checks

BOOST_AUTO_TEST_CASE( test_origin_acceleration) {
    vector<double> pos = {0, 0, 0};

    double r = r_sphere(NFWAcceleration(pos, 1, 10));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));

    r = r_sphere(DiskAcceleration(pos, 1, 0.5));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));
    r = r_sphere(HernquistAcceleration(pos, 1));
    BOOST_TEST(r == 0, tt::tolerance(1e-8));
}


BOOST_DATA_TEST_CASE(test_simple_hernquist, rand_a^bdata::xrange(Npoints), a, index) {
    double r = a;
    vector<double> pos = r * rand_unit();

    double phi_exp = -1. / (2*a);
    double phi = HernquistPotential(pos, a);
    BOOST_TEST(phi==phi_exp, tt::tolerance(1e-6));
    
    double acc_exp = 1 / square(r + a);
    double acc = r_sphere(HernquistAcceleration(pos, a));
    BOOST_TEST(acc==acc_exp, tt::tolerance(1e-6));
}


BOOST_DATA_TEST_CASE(test_simple_nfw, rand_a^rand_c^bdata::xrange(Npoints), R, c, index) {
    double r = R;
    vector<double> pos = r * rand_unit();

    double A = log(1. + c) - c/(1.+c);
    BOOST_TEST(A == NFW_A(c), tt::tolerance(1e-8));

    BOOST_TEST(r_sphere(pos) == R, tt::tolerance(1e-6));
    double phi_exp = -1. / R * log(2) * 1./A;
    double phi = NFWPotential(pos, R, c);
    BOOST_TEST(phi==phi_exp, tt::tolerance(1e-6));
    
    double acc_exp = (1./square(R)) *  1./A * (log(2.) - 1./2.);
    double acc = r_sphere(NFWAcceleration(pos, R, c));
    BOOST_TEST(acc==acc_exp, tt::tolerance(1e-6));
}


BOOST_DATA_TEST_CASE(test_simple_disk, rand_a ^ bdata::xrange(Npoints), a, index) {
    double r = a;
    double b = 0.;
    double z = 0.;
    vector<double> pos =  rand_unit();
    pos[2] = z;
    pos = r * normalize(pos);

    double phi_exp = -1. / (sqrt(2.)*a);
    double phi = DiskPotential(pos, a, b);
    BOOST_TEST(phi==phi_exp, tt::tolerance(1e-6));
    
    double acc_exp = 1. / (square(a) * 2 * sqrt(2));
    double acc_r = r_cyl(DiskAcceleration(pos, a, b));
    double acc_z = z_cyl(DiskAcceleration(pos, a, b));
    BOOST_TEST(acc_r==acc_exp, tt::tolerance(1e-6));
    BOOST_TEST(acc_z==.0, tt::tolerance(1e-6));
}


BOOST_DATA_TEST_CASE(test_simple_disk_z, rand_b ^ bdata::xrange(Npoints), b, index) {
    double a = 0.;
    double z = b;
    fvector pos = {0,0, z};

    double phi_exp = -1. / (sqrt(2.)*b);
    double phi = DiskPotential(pos, a, b);
    BOOST_TEST(phi==phi_exp, tt::tolerance(1e-6));
    
    double acc_exp = 1. / (square(b) * 2 * sqrt(2));
    double acc_r = r_cyl(DiskAcceleration(pos, a, b));
    double acc_z = abs(z_cyl(DiskAcceleration(pos, a, b)));
    BOOST_TEST(acc_r==.0, tt::tolerance(1e-6));
    BOOST_TEST(acc_z==acc_exp, tt::tolerance(1e-6));
}


// sign checks (data)

void test_same_direction(vector<double> a, vector<double>b, double tol=1e-8) {
    vector<double> a_norm = normalize(a);
    vector<double> b_norm = normalize(b);
    double diff = 0;
    for (int i=0; i<3; i++) {
        diff += abs(a_norm[i] - b_norm[i]);
    }
    BOOST_TEST(diff < 1e-8);
}


BOOST_DATA_TEST_CASE(
        test_nfw_acc_sign, 
        bdata::make(vectors()), 
        pos) 
{

    vector<double> acc = NFWAcceleration(pos, 1, 10);
    test_same_direction(acc, -1 * pos);
}


BOOST_DATA_TEST_CASE(
        test_hernquist_acc_sign, 
        bdata::make(vectors()), 
        pos) 
{

    vector<double> acc = HernquistAcceleration(pos, 1);
    test_same_direction(acc, -1 * pos);
}


BOOST_DATA_TEST_CASE(
        test_disk_acc_sign, 
        bdata::make(vectors()), 
        pos) 
{
    vector<double> acc = DiskAcceleration(pos, 4.31, 0.51);
    for (int i=0; i<3; i++) {
        double s = acc[i]*(1*pos)[i];
        BOOST_TEST(s <= 0);
    }
}


BOOST_DATA_TEST_CASE(
        test_nfw_pot_sign, 
        bdata::make(vectors()), 
        pos) 
{

    double pot = NFWPotential(pos, 1, 10);
    BOOST_TEST(pot <= 0);
}


BOOST_DATA_TEST_CASE(
        test_hernquist_pot_sign, 
        bdata::make(vectors()), 
        pos) 
{

    double pot = HernquistPotential(pos, 1);
    BOOST_TEST(pot <= 0);
}

BOOST_DATA_TEST_CASE(
        test_disk_pot_sign, 
        bdata::make(vectors()), 
        pos) 
{

    double pot = DiskPotential(pos, 1, 0.2);
    BOOST_TEST(pot <= 0);
}






// ===========================================
// Derivative tests
// 

struct nfw_params {double R; double c;};

double f_nfw(double r, void * params) {
    struct nfw_params * p = (struct nfw_params *)params;
    double R = (p->R);
    double c = (p->c);
    return NFWPotential(r, R, c);
}


BOOST_DATA_TEST_CASE(
        test_nfw_diff,
        bdata::make(test_radii()) ^ rand_a ^ rand_c,
        r, R, c)
{
    double radius = r + 0.1;
    gsl_function F;
    double result, abserr;
    F.function = &f_nfw;
    struct nfw_params p = {R, c};
    F.params = &p;

    double h = radius * 1e-6;
    gsl_deriv_forward(&F, radius, h, &result, &abserr);
    result = abs(result);

    double actual = radius * NFWScalarAcceleration(radius, p.R, p.c);
    BOOST_TEST(result == actual, tt::tolerance(1e-4));
}


double f_hernquist(double r, void * params) {
    double a = (*(double*) params);
    return HernquistPotential(r, a);
}


BOOST_DATA_TEST_CASE(
        test_hernquist_ddiff,
        bdata::make(test_radii()) ^ rand_a,
        r, a)
{
    double radius = r + 0.1;
    double A = a;
    gsl_function F;
    double result, abserr;
    F.function = &f_hernquist;
    F.params = &A;

    double h = radius * 1e-6;
    gsl_deriv_forward(&F, radius, h, &result, &abserr);
    result = abs(result);

    double actual = radius * HernquistScalarAcceleration(radius, a);
    BOOST_TEST(result == actual, tt::tolerance(1e-4));
}



struct disk_params {double R; double a; double b; double z;};

double f_disk_z(double z, void * params) {
    struct disk_params * p = (struct disk_params *)params;
    double a = (p->a);
    double b = (p->b);
    double R = (p->R);
    return DiskPotential(R, z, a, b);
}


double f_disk_R(double R, void * params) {
    struct disk_params * p = (struct disk_params *)params;
    double a = (p->a);
    double b = (p->b);
    double z = (p->z);
    return DiskPotential(R, z, a, b);
}


BOOST_DATA_TEST_CASE(
        test_disk_diff_z,
        bdata::make(vectors()) ^ rand_a ^ rand_b,
        pos, a, b)
{
    double r = r_sphere(pos);
    fvector pos1 = (r+1) * normalize(pos);

    double R = r_cyl(pos1);
    double z = z_cyl(pos1);
    struct disk_params p = {R, a, b, z};

    gsl_function F;
    F.params = &p;

    double result_R, abserr;
    double h = R * 1e-6;
    F.function = &f_disk_R;
    gsl_deriv_forward(&F, R, h, &result_R, &abserr);

    double result_z;
    h = z * 1e-6;
    F.function = &f_disk_z;
    gsl_deriv_forward(&F, z, h, &result_z, &abserr);


    fvector acc =  DiskAcceleration(pos1, a, b);
    double acc_R = r_cyl(acc);
    double acc_z = z_cyl(acc);
    BOOST_TEST(acc_z == -result_z, tt::tolerance(1e-4));
    BOOST_TEST(acc_R == result_R, tt::tolerance(1e-4));
}

