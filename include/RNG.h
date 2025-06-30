/*
 *  RNG.h
 *
 *  DEBUG VERSION
 *  Header file containing the declaration of the RNG class including
 *  the full definitions of the trivial member functions ("getters").
 *  The full definition of the class is found in the associated C++ file
 *  "RNG.cpp"
 *
 *  James Bradford
 *  University of Sheffield
 *  Apr 2017
 */

#ifndef RNG_H
#define RNG_H
#include <array>

class RNG
{
private:
    std::random_device m_seed_generator;
    std::uint32_t m_random_seed;
public:
    // Default Constructor
    RNG();

    std::mt19937_64 m_mersenne;
    std::array<std::mt19937_64, 16> m_PRNGs;

    std::uniform_real_distribution<double> m_angledist;
    std::normal_distribution<double> m_posdist;
    std::uniform_real_distribution<double> m_probDist;

    std::uint32_t getSeed() const { return m_random_seed; }
    auto getRN() const { return m_mersenne; }

    void setSeed(std::uint32_t seed);

};

#endif
