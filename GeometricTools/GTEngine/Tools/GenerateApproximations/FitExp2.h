// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngine.h>
using namespace gte;

class FitExp2
{
public:
    FitExp2();

    template <int Degree>
    void Generate(std::vector<double>& poly, double& error) const;

private:
    // f(x) = 2^x for x in [0,1].
    // g(x) = 2^x - p(x)
    double G(double x, int const degree, double const* p) const;

    // g'(x) = 2^x*log(2) - p'(x)
    double GDer(double x, int const degree, double const* p) const;

    int SignG(double x, int const degree, double const* p) const;
    int SignGDer(double x, int const degree, double const* p) const;
};


FitExp2::FitExp2()
{
}

double FitExp2::G(double x, int const degree, double const* p) const
{
    int i = degree;
    double result = p[i];
    while (--i >= 0)
    {
        result = x*result + p[i];
    }
    result = pow(2.0, x) - result;
    return result;
}

double FitExp2::GDer(double x, int const degree, double const* p) const
{
    int i = degree;
    double result = i*p[i];
    while (--i >= 1)
    {
        result = x*result + i*p[i];
    }
    result = log(2.0)*pow(2.0, x) - result;
    return result;
}

int FitExp2::SignG(double x, int const degree, double const* p) const
{
    double g = G(x, degree, p);
    return (g > 0.0 ? 1 : (g < 0.0 ? -1 : 0));
}

int FitExp2::SignGDer(double x, int const degree, double const* p) const
{
    double gder = GDer(x, degree, p);
    return (gder > 0.0 ? 1 : (gder < 0.0 ? -1 : 0));
}

template <int Degree>
void FitExp2::Generate(std::vector<double>& poly, double& error) const
{
    Vector<Degree + 1, double> g0root;
    g0root[0] = 0.0;
    for (int i = 1; i < Degree; ++i)
    {
        g0root[i] = static_cast<double>(i) / static_cast<double>(Degree);
    }
    g0root[Degree] = 1.0;

    // Initial guess at the polynomial, gives you oscillation.
    Matrix<Degree + 1, Degree + 1, double> A;
    Vector<Degree + 1, double> B;
    for (int r = 0; r <= Degree; ++r)
    {
        A(r, 0) = 1.0;
        for (int c = 1; c <= Degree; ++c)
        {
            A(r, c) = g0root[r] * A(r, c - 1);
        }
        B[r] = pow(2.0, g0root[r]);
    }
    Vector<Degree + 1, double> p = Inverse(A) * B;
    double absError[Degree], e;

    for (int iter = 0; iter < 16; ++iter)
    {
        // Find local extremum on [g0root[i],g0root[i+1]].
        Vector<Degree, double> g1root;
        int k = 0;
        for (int i = 0; i < Degree; ++i)
        {
            double x0 = g0root[i];
            double x1 = g0root[i + 1];
            int s0 = SignGDer(x0, Degree, &p[0]);
            int s1 = SignGDer(x1, Degree, &p[0]);
            LogAssert(s0*s1 < 0, "Unexpected condition.");
            int j;
            for (j = 0; j < 1024; ++j)
            {
                double xmid = 0.5*(x0 + x1);
                if (x0 != xmid && x1 != xmid)
                {
                    int smid = SignGDer(xmid, Degree, &p[0]);
                    if (smid == s0)
                    {
                        x0 = xmid;
                    }
                    else if (smid == s1)
                    {
                        x1 = xmid;
                    }
                    else
                    {
                        g1root[k++] = xmid;
                        break;
                    }
                }
                else
                {
                    g1root[k++] = xmid;
                    break;
                }

            }
            LogAssert(j < 1024, "Max j-iterations not large enough.");
        }

        for (int i = 0; i < Degree; ++i)
        {
            absError[i] = std::abs(G(g1root[i], Degree, &p[0]));
        }

        // Solve g(x[i]) + (-1)^{i}*e = f(x[i]) for e and p[j].
        double sign = 1.0;
        for (int r = 0; r < Degree; ++r, sign = -sign)
        {
            A(r, 0) = g1root[r];
            for (int c = 1; c < Degree; ++c)
            {
                A(r, c) = g1root[r] * A(r, c - 1);
            }
            A(r, Degree) = sign;
            B[r] = pow(2.0, g1root[r]) - 1.0;
        }
        for (int c = 0; c < Degree; ++c)
        {
            A(Degree, c) = 1.0;
        }
        A(Degree, Degree) = 0.0;
        B[Degree] = 2.0 - 1.0;

        Vector<Degree + 1, double> X = Inverse(A)*B;
        for (int i = 0; i < Degree; ++i)
        {
            p[i + 1] = X[i];
        }
        e = X[Degree];

        k = 1;
        for (int i = 0; i < Degree - 1; ++i)
        {
            double x0 = g1root[i];
            double x1 = g1root[i + 1];
            int s0 = SignG(x0, Degree, &p[0]);
            int s1 = SignG(x1, Degree, &p[0]);
            LogAssert(s0*s1 < 0, "Unexpected condition.");
            int j;
            for (j = 0; j < 1024; ++j)
            {
                double xmid = 0.5*(x0 + x1);
                if (x0 != xmid && x1 != xmid)
                {
                    int smid = SignG(xmid, Degree, &p[0]);
                    if (smid == s0)
                    {
                        x0 = xmid;
                    }
                    else if (smid == s1)
                    {
                        x1 = xmid;
                    }
                    else
                    {
                        g0root[k++] = xmid;
                        break;
                    }
                }
                else
                {
                    g0root[k++] = (x0 == xmid ? x0 : x1);
                    break;
                }

            }
            LogAssert(j < 1024, "Max j-iterations not large enough.");
        }
    }

    poly.resize(Degree + 1);
    for (int i = 0; i <= Degree; ++i)
    {
        poly[i] = p[i];
    }
    error = e;
}

