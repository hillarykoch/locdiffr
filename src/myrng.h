#ifndef _MYRNG_H
#define _MYRNG_H

#include <random>
using namespace std;

class Uniform {
    public:
        mt19937 _generator;
        uniform_real_distribution<double> _my_unif;
        Uniform(unsigned int seed) {
            _generator = mt19937(seed);
            _my_unif = uniform_real_distribution<>(0.0, 1.0);
        }

        double draw(void) {
            return _my_unif(_generator);
        }
};

#endif
