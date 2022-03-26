#ifndef WEIGHTS_HPP
#define WEIGHTS_HPP

#include "cdeepso_params.hpp"

#include <wup/wup.hpp>
#include <vector>

using namespace wup;
using namespace std;

class Weight
{
public:

    Random * generator;

    Precision pInertia; // 1
    Precision pMemory; // 2
    Precision pCooperation; // 3
    Precision pPerturbation; // 4
    Precision dThreshold; // 5
    Precision dVelocity; // 6

public:

    Weight() {}

    void
    init(Random & generator,
         Precision const maxVelocity)
    {
        this->generator = & generator;
        pInertia = generator.uniformDouble();
        pMemory = generator.uniformDouble();
        pCooperation = generator.uniformDouble();
        pPerturbation = generator.uniformDouble();
        dThreshold = generator.uniformDouble();
        dVelocity = generator.uniformDouble() * maxVelocity;
    }

    void
    copyWithNoise(Weight const & s,
                  Precision const mutationRate,
                  Precision const maxVelocity)
    {
        pInertia = copyWithNoise(s.pInertia, mutationRate, 1.0);
        pMemory = copyWithNoise(s.pMemory, mutationRate, 1.0);
        pCooperation = copyWithNoise(s.pCooperation, mutationRate, 1.0);
        pPerturbation = copyWithNoise(s.pPerturbation, mutationRate, 1.0);
        dThreshold = copyWithNoise(s.dThreshold, mutationRate, 1.0);
        dVelocity = copyWithNoise(s.dVelocity, mutationRate, maxVelocity);
//        print("weight copied :", pInertia, pMemory, pCooperation, pPerturbation, dThreshold, dVelocity);
    }

    Precision
    copyWithNoise(Precision const s,
                  Precision const mutationRate,
                  Precision const max)
    {
        // Method 1
        Precision v = s + generator->normalDouble() * mutationRate;
        if (v < 0.0) return 0.0;
        if (v > max) return max;
        return v;

        // Method 2
//        Precision v = s + generator->gaussianNoise() * p.mutationRate;
//        return std::fmod(v,max);

        // Method 3
//        Precision v = s + generator->gaussianNoise() * p.mutationRate;
//        return abs((fmod(v, 2.0) - 0.5)) * max;

        // Method 4
//        return generator->uniformNoise() * max;
    }

};

#endif // WEIGHTS_HPP
