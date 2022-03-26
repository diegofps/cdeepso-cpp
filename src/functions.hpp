#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "cdeepso_params.hpp"
#include "population.hpp"
#include "utils.hpp"


//////////////////////////////////////////////////////////////////////////////////////////
// Particle functions
//////////////////////////////////////////////////////////////////////////////////////////

Precision
rosenbrock(const Precision * const x, const int len)
{
    Precision fit = 0.0;

    for (int i=0;i!=len-1;++i)
    {
        Precision term1 = x[i+1] - x[i]*x[i];
        Precision term2 = 1 - x[i];
        fit += 100 * term1*term1 + term2*term2;
    }

    return fit;
}

Precision
griewank(const Precision * const x, const int len)
{
    Precision term1 = 0;
    Precision term2 = 1;

    for (int i=0;i!=len;++i)
    {
        term1 += x[i]*x[i];
        term2 = term2 * cos( x[i] / sqrt(i+1) );
    }

    return 1 + term1 / 4000 - term2;
}

Precision
rastrigin(const Precision * const x, const int len)
{
    double fit = 0.0;

    for (int i=0;i!=len;++i)
        fit += x[i]*x[i] - 10 * cos(2 * M_PI * x[i]);

    return 10 * len + fit;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Population functions
//////////////////////////////////////////////////////////////////////////////////////////

void
rosenbrock(Particles & particles,
           Refreshes & refresh,
           Fitness & fitness)
{
    for (uint i=0;i!=particles.numRows();++i)
        if (refresh[i])
            fitness[i] = rosenbrock(&particles(i,0), particles.numCols());
}

void
griewank(Particles & particles,
         Refreshes & refresh,
         Fitness & fitness)
{
    for (uint i=0;i!=particles.numRows();++i)
        if (refresh[i])
            fitness[i] = griewank(&particles(i,0), particles.numCols());
}

void
rastrigin(Particles & particles,
          Refreshes & refresh,
          Fitness & fitness)
{
    for (uint i=0;i!=particles.numRows();++i)
        if (refresh[i])
            fitness[i] = rastrigin(&particles(i,0), particles.numCols());
}


//////////////////////////////////////////////////////////////////////////////////////////
// Other functions
//////////////////////////////////////////////////////////////////////////////////////////


Precision
eval1(const Precision * const particle, const int len)
{
    return wup::arr::max(particle, len);
}


Precision
eval2(const Precision * const particle, const int len)
{
    Precision sum = 0.0;
    for (int i=0;i!=len;++i)
        sum += particle[i];
    return sum;
}


class Rosenbrock
{
public:
    vector<Precision> term1;
    vector<Precision> term2;
public:

    Rosenbrock(CDEEPSOParams & p) :
        term1(p.dims-1), term2(p.dims-1) { }

    Precision operator()(const Precision * const x, const int len)
    {
        for (int i=0;i!=len-1;++i)
        {
            term1[i] = x[i+1] - x[i]*x[i];
            term2[i] = 1.0 - x[i];
        }

        arraySquaredCollapse(term1);
        arraySquaredCollapse(term2);

        return 100 * arraySumCollapse(term1) + arraySumCollapse(term2);
    }
};


class Rastrigin
{
public:
    vector<Precision> tmp1;
    vector<Precision> tmp2;

public:
    Rastrigin(CDEEPSOParams & p) :
        tmp1(p.dims), tmp2(p.dims) { }

    Precision
    operator()(const Precision * const x, const int len)
    {
        for (int i=0;i!=len;++i)
        {
            tmp1[i] = x[i];
            tmp2[i] = cos(2 * M_PI * x[i]);
        }

        arraySquaredCollapse(tmp1);

        // x[i]*x[i] - 10 * cos(2 * M_PI * x[i]);

        return arraySumCollapse(tmp1) - 10 * arraySumCollapse(tmp2);
    }
};


#endif // FUNCTIONS_HPP
