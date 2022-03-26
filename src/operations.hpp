#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "population.hpp"

namespace ops
{

inline void
initPopulation(Population & current,
               Random & generator,
               vector<double> const & xMin,
               vector<double> const & xMax,
               vector<double> const & vMin,
               vector<double> const & vMax,
               Precision const maxVelocity)
{
    for (uint i=0;i!=current.size();++i)
        current.weights[i].init(generator, maxVelocity);

    for (uint i=0;i!=current.size();++i)
    {
        for (uint j=0;j!=current.dims();++j)
        {
            current.particles(i,j) = xMin[j] + (xMax[j] - xMin[j]) * generator.uniformDouble();
            current.velocity(i,j) = vMin[j] + (vMax[j] - vMin[j]) * generator.uniformDouble();
        }
    }
}

inline void
initBests(Population & current,
          Fitness & fitness,
          Population & myBest,
          Fitness & myBestFitness,
          vector<Precision> & gBest,
          Precision & gBestFit)
{
    myBest.cloneFrom(current);
    myBestFitness = fitness;
    const int srcId = arr::indexOfMin(fitness);
    current.particles.exportRow(srcId, gBest);
    gBestFit = fitness[srcId];
}

inline void
initLimits(int const dims,
           int const minValue,
           int const maxValue,
           vector<double> & xMin,
           vector<double> & xMax,
           vector<double> & vMin,
           vector<double> & vMax)
{
    for (int i=0;i!=dims;++i)
    {
        xMin[i] = minValue;
        xMax[i] = maxValue;
        vMin[i] = xMin[i] - xMax[i];
        vMax[i] = -vMin[i];
    }
}


inline void
computeNewWeights(const Population & src,
                  Population & dst,
                  Precision const mutationRate,
                  Precision const maxVelocity)
{
    for (uint i=0;i!=src.size();++i)
        dst.weights[i].copyWithNoise(src.weights[i], mutationRate, maxVelocity);
}

inline void
computeNewVel(Population & pop,
              Random & generator,
              Population const & myBest,
              vector<Precision> const & gBest,
              vector<double> const & vMin,
              vector<double> const & vMax,
              Precision const communicationProbability)
{
    for (uint i=0;i!=pop.size();++i)
    {
        const Weight & weight = pop.weights[i];
        const Precision * pos = & pop.particles(i,0);
        const Precision * vel = & pop.velocity(i,0);
        const Precision * mbp = & myBest.particles(i,0);

        Precision * d = & pop.velocity(i,0);
        const Precision noise = 1.0 + weight.pPerturbation * generator.normalDouble();

        for (uint k=0;k!=pop.dims();++k)
        {
            const double it = weight.pInertia * vel[k];
            const double mt = weight.pMemory * (mbp[k] - pos[k]);

            const Precision ct = generator.unfairCoin(communicationProbability)
                    ? weight.pCooperation * (gBest[k] * noise - pos[k])
                    : 0.0;

            const double tmp = it + mt + ct;

            d[k] = tmp > vMax[k] ? vMax[k] : tmp < vMin[k] ? vMin[k] : tmp;
        }
    }
}

inline void
computeNewPos(Population & pop)
{
    const Precision * const posEnd = pop.particles.end();
    const Precision * posPtr = pop.particles.begin();
    const Precision * velPtr = pop.velocity.begin();
    Precision * newPosPtr = pop.particles.begin();

    while(posPtr != posEnd)
    {
        *newPosPtr = *posPtr + *velPtr;

        ++newPosPtr;
        ++posPtr;
        ++velPtr;
    }
}

inline void
enforceLimits(Population & pop,
              vector<double> & xMin,
              vector<double> & xMax,
              vector<double> & vMin,
              vector<double> & vMax)
{
    for (uint i=0;i!=pop.size();++i)
    {
        for (uint j=0;j!=pop.dims();++j)
        {
            if (pop.particles(i,j) < xMin[j])
            {
                pop.particles(i,j) = xMin[j];
                if (pop.velocity(i,j) < 0)
                    pop.velocity(i,j) = -pop.velocity(i,j);
            }

            else if (pop.particles(i,j) > xMax[j])
            {
                pop.particles(i,j) = xMax[j];
                if (pop.velocity(i,j) > 0)
                    pop.velocity(i,j) = -pop.velocity(i,j);
            }

            if (pop.velocity(i,j) < vMin[j])
                pop.velocity(i,j) = vMin[j];

            else if (pop.velocity(i,j) > vMax[j])
                pop.velocity(i,j) = vMax[j];
        }
    }
}

inline void
mergePopulations(Population const & src,
                 Population & dst,
                 Fitness & srcFitness,
                 Fitness & dstFitness)
{
    for (uint i=0;i!=src.size();++i)
    {
        if (srcFitness[i] < dstFitness[i])
        {
            dst.particles.importRow(src.particles, i, i);
            dst.velocity.importRow(src.velocity, i, i);
            dst.weights[i] = src.weights[i];
//            dstFitness[i] = srcFitness[i];
        }
    }
}

inline void
updateGBest(Population const & pop,
            Fitness const & popFitness,
            Population & memGBest,
            Fitness & memGBestFitness,
            int & memGBestIndex,
            vector<Precision> & gBest,
            Precision & gBestFit)
{
    const int srcId = arr::indexOfMin(popFitness);

    if (popFitness[srcId] < gBestFit)
    {
        pop.particles.exportRow(srcId, gBest);
        gBestFit = popFitness[srcId];

        const int dstId = memGBestIndex == int(memGBest.size())
                ? arr::indexOfMax(memGBestFitness)
                : memGBestIndex++;

        memGBest.particles.importRow(pop.particles, srcId, dstId);
        memGBest.velocity.importRow(pop.velocity, srcId, dstId);
        memGBest.weights[dstId] = pop.weights[srcId];
        memGBestFitness[dstId] = popFitness[srcId];
    }
}

inline void
updateMyBestPos(Population const & pop,
                Fitness const & popFitness,
                Population & myBest,
                Fitness & myBestFitness)
{
    for (uint i=0;i!=pop.size();++i)
    {
        if (popFitness[i] < myBestFitness[i])
        {
            myBest.particles.importRow(pop.particles, i, i);
            myBest.velocity.importRow(pop.velocity, i, i);
            myBest.weights[i] = pop.weights[i];
            myBestFitness[i] = popFitness[i];
        }
    }
}

inline void
updateCandidates(int const k,
                 Population const & pop,
                 Fitness const & popFitness,
                 Fitness & memGBestFitness,
                 vector<int> & candidates,
                 int const memGBestIndex,
                 CDEEPSOParams::MemStrategy const memStrategy)
{
    const double particleFit = popFitness[k];
    candidates.clear();

    if (memStrategy & CDEEPSOParams::MemStrategy::MEM)
        for (int i=0;i!=memGBestIndex;++i)
            if (memGBestFitness[i] < particleFit)
                candidates.push_back(-i);

    if (memStrategy & CDEEPSOParams::MemStrategy::POS)
        for (uint i=0;i!=pop.size();++i)
            if (popFitness[i] < particleFit)
                candidates.push_back(i+1);
}

inline void
heuristicRand(Population const & src,
              Fitness const & srcFitness,
              Population & dst,
              Population & myBest,
              Population & memGBest,
              Fitness & memGBestFitness,
              int const memGBestIndex,
              CDEEPSOParams::MemStrategy const memStrategy,
              vector<int> & candidates,
              vector<bool> & dstRefresh,
              Random & generator)
{
    for (uint i=0;i!=src.size();++i)
    {
        updateCandidates(i, src, srcFitness, memGBestFitness, candidates, memGBestIndex, memStrategy);

        if (candidates.size() >= 3)
        {
            generator.shuffle(candidates);

            Precision const * const mgb1 = candidates[0] > 0
                    ? & src.particles(candidates[0]-1,0)
                    : & memGBest.particles(-candidates[0],0);

            Precision const * const mgb2 = candidates[1] > 0
                    ? & src.particles(candidates[1]-1,0)
                    : & memGBest.particles(-candidates[1],0);

            Precision const * const mgb3 = candidates[2] > 0
                    ? & src.particles(candidates[2]-1,0)
                    : & memGBest.particles(-candidates[2],0);

            Precision const * const mbp = & myBest.particles(i,0);
            Weight const & w = src.weights[i];
            Precision * const d = & dst.particles(i,0);
            dstRefresh[i] = true;
            dst.weights[i] = w;

            for (uint j=0;j!=src.dims();++j)
                d[j] = mgb1[j] + w.dVelocity * (mgb2[j] - mgb3[j]);

            uint const tmpIndexD = generator.uniformInt(src.dims());

            for (uint j=0;j!=src.dims();++j)
                if (generator.unfairCoin(w.dThreshold) || j == tmpIndexD)
                    d[j] = mbp[j];
        }
        else
        {
            dst.particles.importRow(src.particles, i, i);
            dst.velocity.importRow(src.velocity, i, i);
            dst.weights[i] = src.weights[i];
        }
    }
}

void
heuristicBest(const Population & src,
              Fitness const & srcFitness,
              Population & dst,
              vector<Precision> & gBest,
              Population & memGBest,
              Fitness & memGBestFitness,
              int const memGBestIndex,
              CDEEPSOParams::MemStrategy const memStrategy,
              vector<int> & candidates,
              vector<bool> & dstRefresh,
              Random & generator)
{
    for (uint i=0;i!=src.size();++i)
    {
        updateCandidates(i, src, srcFitness, memGBestFitness, candidates, memGBestIndex, memStrategy);

        if (candidates.size() >= 2)
        {
            generator.shuffle(candidates);

            Precision const * const mgb1 = candidates[0] > 0
                    ? & src.particles(candidates[0]-1,0)
                    : & memGBest.particles(-candidates[0],0);

            Precision const * const mgb2 = candidates[1] > 0
                    ? & src.particles(candidates[1]-1,0)
                    : & memGBest.particles(-candidates[1],0);

            Weight const & w = src.weights[i];
            Precision * const d = & dst.particles(i,0);
            dstRefresh[i] = true;
            dst.weights[i] = w;

            for (uint j=0;j!=src.dims();++j)
                d[j] = gBest[j] + w.dVelocity * (mgb1[j] - mgb2[j]);

            uint const tmpIndexD = generator.uniformInt(src.dims());

            for (uint j=0;j!=src.dims();++j)
                if (generator.unfairCoin(w.dThreshold) || j == tmpIndexD)
                    d[j] = gBest[j];
        }
        else
        {
            dst.particles.importRow(src.particles, i, i);
            dst.velocity.importRow(src.velocity, i, i);
            dst.weights[i] = src.weights[i];
        }
    }
}

}

#endif // OPERATIONS_HPP
