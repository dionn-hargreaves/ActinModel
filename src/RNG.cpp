/*
 *  RNG.cpp
 *
 *  C++ file containing the definition of the Monte Carlo class.
 *  Here is the home of the RNG, and associated distributions
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "RNG.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Apr 2017
 */

#include "configHeader.h"
#include "RNG.h"

 /*
  * RANDOM NUMBER GENERATER
  * The generator itself is a static variable, as we only want to seed it once
  * for every run!
  * The distributions are always static, as they will be the same for all objects
  * therefore it makes sense to save memory.
  * Remember static members are created when the program starts and are
  * associated with the class not any object
  * random_device is a non-deterministic uniform rng using a hardware entropy
  * source.
  * mt19937_64 generates 64-bit unsigned integers using the output from
  * random_device as a seed.
  */


RNG::RNG() //Constructor 0
{
    // Define distributions that will be used to post-process the output from the
    // rng we defined above.
    m_angledist = std::uniform_real_distribution<double> (0, (2*M_PI));
    m_posdist = std::normal_distribution<double> (0,1E-6);
    m_probDist = std::uniform_real_distribution<double> (0, 1);
}

void RNG::setSeed(std::uint32_t seed)
{
    m_random_seed = seed;
    m_mersenne = std::mt19937_64 (m_random_seed);
    m_PRNGs = { std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne),
              std::mt19937_64 (m_mersenne)  };
}
