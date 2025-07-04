// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteConstants.h>
#include <cmath>

// Minimax polynomial approximations to 2^x.  The polynomial p(x) of
// degree D minimizes the quantity maximum{|2^x - p(x)| : x in [0,1]}
// over all polynomials of degree D.

namespace gte
{

template <typename Real>
class Exp2Estimate
{
public:
    // The input constraint is x in [0,1].  For example,
    //   float x; // in [0,1]
    //   float result = Exp2Estimate<float>::Degree<3>(x);
    template <int D>
    inline static Real Degree(Real x);

    // The input x can be any real number.  Range reduction is used to
    // generate a value y in [0,1], call Degree(y), and combine the output
    // with the proper exponent to obtain the approximation.  For example,
    //   float x;  // x >= 0
    //   float result = Exp2Estimate<float>::DegreeRR<3>(x);
    template <int D>
    inline static Real DegreeRR(Real x);

private:
    // Metaprogramming and private implementation to allow specialization of
    // a template member function.
    template <int D> struct degree {};
    inline static Real Evaluate(degree<1>, Real t);
    inline static Real Evaluate(degree<2>, Real t);
    inline static Real Evaluate(degree<3>, Real t);
    inline static Real Evaluate(degree<4>, Real t);
    inline static Real Evaluate(degree<5>, Real t);
    inline static Real Evaluate(degree<6>, Real t);
    inline static Real Evaluate(degree<7>, Real t);
};


template <typename Real>
template <int D>
inline Real Exp2Estimate<Real>::Degree(Real x)
{
    return Evaluate(degree<D>(), x);
}

template <typename Real>
template <int D>
inline Real Exp2Estimate<Real>::DegreeRR(Real x)
{
    Real p = floor(x);
    Real y = x - p;
    Real poly = Degree<D>(y);
    Real result = ldexp(poly, (int)p);
    return result;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<1>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG1_C1;
    poly = (Real)GTE_C_EXP2_DEG1_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<2>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG2_C2;
    poly = (Real)GTE_C_EXP2_DEG2_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG2_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<3>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG3_C3;
    poly = (Real)GTE_C_EXP2_DEG3_C2 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG3_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG3_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<4>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG4_C4;
    poly = (Real)GTE_C_EXP2_DEG4_C3 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG4_C2 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG4_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG4_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<5>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG5_C5;
    poly = (Real)GTE_C_EXP2_DEG5_C4 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG5_C3 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG5_C2 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG5_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG5_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<6>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG6_C6;
    poly = (Real)GTE_C_EXP2_DEG6_C5 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG6_C4 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG6_C3 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG6_C2 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG6_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG6_C0 + poly * t;
    return poly;
}

template <typename Real>
inline Real Exp2Estimate<Real>::Evaluate(degree<7>, Real t)
{
    Real poly;
    poly = (Real)GTE_C_EXP2_DEG7_C7;
    poly = (Real)GTE_C_EXP2_DEG7_C6 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C5 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C4 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C3 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C2 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C1 + poly * t;
    poly = (Real)GTE_C_EXP2_DEG7_C0 + poly * t;
    return poly;
}


}
