#ifndef POPULATION_HPP
#define POPULATION_HPP

#include "weight.hpp"

typedef Bundle<Precision> Particles;
typedef Bundle<Precision> Velocities;
typedef vector<Weight> Weights;
typedef vector<Precision> Fitness;
typedef vector<bool> Refreshes;

class Population
{
public:

    Particles particles;
    Velocities velocity;
    Weights weights;

public:

    Population(const int popSize, const int dims) :
        particles(popSize, dims, 0),
        velocity(popSize, dims, 0),
        weights(popSize)
    {

    }

    void
    cloneFrom(const Population & other)
    {
        for (size_t i=0;i!=particles.numRows();++i)
            weights[i] = other.weights[i];

        for (size_t i=0;i!=particles.numRows();++i)
            for (size_t j=0;j!=particles.numCols();++j)
            {
                particles(i,j) = other.particles(i,j);
                velocity(i,j) = other.velocity(i,j);
            }

//        copy(other.pos.begin(), other.pos.end(), pos.begin());
//        copy(other.vel.begin(), other.vel.end(), vel.begin());
//        copy(other.weights.begin(), other.weights.end(), weights.begin());
//        copy(other.fitness.begin(), other.fitness.end(), fitness.begin());
    }

    uint
    size() const
    {
        return particles.numRows();
    }

    uint
    dims() const
    {
        return particles.numCols();
    }

};

#endif // POPULATION_HPP
