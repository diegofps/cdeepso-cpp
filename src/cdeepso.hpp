#ifndef CDEEPSO_HPP
#define CDEEPSO_HPP

#include "cdeepso_params.hpp"
#include "operations.hpp"
#include "population.hpp"
#include "weight.hpp"

class CDEEPSO
{
public:

    CDEEPSOParams & p;

    vector<double> xMin;
    vector<double> xMax;
    vector<double> vMin;
    vector<double> vMax;

    Population pop1;
    Population myBest;
    Population pop2;
    Population memGBest;
    Fitness myBestFitness;
    Fitness memGBestFitness;

    Precision gBestFit;
    vector<Precision> gBest;

    int memGBestIndex;
    vector<int> candidates;
    int fitEval;

    Random generator;

    typedef std::function<void(int const generation, CDEEPSO&)> LoopListener;
    LoopListener onLoopListener;

public:

    CDEEPSO(CDEEPSOParams & p) :
        p(p),

        xMin(p.dims),
        xMax(p.dims),
        vMin(p.dims),
        vMax(p.dims),

        pop1(p.popSize, p.dims),
        myBest(p.popSize, p.dims),
        pop2(p.popSize, p.dims),
        memGBest(p.popSize, p.dims),

        myBestFitness(p.popSize),
        memGBestFitness(p.popSize),

        gBestFit(-1.0),
        gBest(p.dims),

        memGBestIndex(0),
        fitEval(0)

    {
        ops::initLimits(p.dims, p.xMin, p.xMax, xMin, xMax, vMin, vMax);
        candidates.reserve(p.popSize + p.memGBestSize);
    }

    void
    setOnLoopListener(LoopListener onLoopListener)
    {
        this->onLoopListener = onLoopListener;
    }

    void
    initPopulationInPop1()
    {
        ops::initPopulation(pop1, generator, xMin, xMax, vMin, vMax, p.maxVelocity);
    }

    void
    initBestsFromPop1(Fitness & pop1Fitness)
    {
        ops::initBests(pop1, pop1Fitness, myBest, myBestFitness, gBest, gBestFit);
    }

    void
    createPop2FromHeuristic(Fitness & pop1Fitness,
                            Refreshes & pop2Refresh)
    {
        if (p.deType == CDEEPSOParams::DEType::RAND)
            ops::heuristicRand(pop1, pop1Fitness, pop2, myBest, memGBest, memGBestFitness, memGBestIndex, p.memStrategy, candidates, pop2Refresh, generator);

        else if (p.deType == CDEEPSOParams::DEType::BEST)
            ops::heuristicBest(pop1, pop1Fitness, pop2, gBest, memGBest, memGBestFitness, memGBestIndex, p.memStrategy, candidates, pop2Refresh, generator);

        else
            error("Unknown deType");

        ops::enforceLimits(pop2, xMin, xMax, vMin, vMax);
    }

    void
    createPop2FromMutatedWeight()
    {
        pop2.cloneFrom(pop1);
        ops::computeNewWeights(pop1, pop2, p.mutationRate, p.maxVelocity);
        ops::computeNewVel(pop2, generator, myBest, gBest, vMin, vMax, p.communicationProbability);
        ops::computeNewPos(pop2);
        ops::enforceLimits(pop2, xMin, xMax, vMin, vMax);
    }

    void
    createPop1FromVelocity()
    {
        ops::computeNewVel(pop1, generator, myBest, gBest, vMin, vMax, p.communicationProbability);
        ops::computeNewPos(pop1);
        ops::enforceLimits(pop1, xMin, xMax, vMin, vMax);
    }

    void
    mergeIntoPop1(Fitness & pop1Fitness,
                  Fitness & pop2Fitness)
    {
        ops::mergePopulations(pop2, pop1, pop2Fitness, pop1Fitness);
        ops::updateMyBestPos(pop1, pop1Fitness, myBest, myBestFitness);
        ops::updateGBest(pop1, pop1Fitness, memGBest, memGBestFitness, memGBestIndex, gBest, gBestFit);
    }

    template <typename EVAL>
    void
    computeFitness(Population & pop,
                   Refreshes & refresh,
                   Fitness & fitness,
                   EVAL eval)
    {
        eval(pop.particles, refresh, fitness);

//        const int oldFitEval = fitEval;

        for (uint i=0;i!=pop.size();++i)
        {
            if (refresh[i])
            {
                fitEval += 1;
                refresh[i] = false;
            }
        }

//        print("Added evals:", fitEval - oldFitEval);
    }

    void
    clearRefresh(Refreshes & refresh,
                 bool const value)
    {
        for (size_t i=0;i!=refresh.size();++i)
            refresh[i] = value;
    }

    template <typename EVAL>
    void
    optimize(EVAL eval, bool initPop=true)
    {
        Fitness pop1Fitness(pop1.size());
        Fitness pop2Fitness(pop2.size());

        Refreshes pop1Refresh(pop1.size());
        Refreshes pop2Refresh(pop2.size());

        int i;

        if (initPop)
            initPopulationInPop1();

        clearRefresh(pop1Refresh, true);
        computeFitness(pop1, pop1Refresh, pop1Fitness, eval);
        initBestsFromPop1(pop1Fitness);

        for (i=0;i!=p.maxGen && fitEval<=p.maxFitEval;++i)
        {
            pop2Fitness = pop1Fitness;
            clearRefresh(pop2Refresh, false);
            createPop2FromHeuristic(pop1Fitness, pop2Refresh);
            computeFitness(pop2, pop2Refresh, pop2Fitness, eval);
            mergeIntoPop1(pop1Fitness, pop2Fitness);

            createPop2FromMutatedWeight();
            clearRefresh(pop2Refresh, true);
            computeFitness(pop2, pop2Refresh, pop2Fitness, eval);

            createPop1FromVelocity();
            clearRefresh(pop1Refresh, true);
            computeFitness(pop1, pop1Refresh, pop1Fitness, eval);

            mergeIntoPop1(pop1Fitness, pop2Fitness);

            if (p.printConvergenceResults != 0 && i % p.printConvergenceResults == 0)
                printn(BLUE, "Gen: ", i, ", Best Fit: ", std::scientific, gBestFit, std::defaultfloat, ", fitEvals:" , fitEval, "/", p.maxFitEval, "\n", NORMAL);

            if (onLoopListener)
                onLoopListener(i, *this);
        }

        printn(YELLOW, "Optimization has ended, Generations: ", i, ", Best Fit: ", std::scientific, gBestFit, std::defaultfloat, ", Fit Evals:" , fitEval, "/", p.maxFitEval, "\n", NORMAL);

    }

public:

    void
    showPop1()
    {
        showPopulation(pop1, p.popSize, "POP1");
    }

    void
    showPop2()
    {
        showPopulation(pop2, p.popSize, "POP2");
    }

    void
    showMemory()
    {
        showPopulation(memGBest, memGBestIndex, "MEMORY");
    }

    void
    showMyBest()
    {
        showPopulation(myBest, p.popSize, "MY_BEST");
    }

    void
    showGBest()
    {
        print(BLUE, "GBEST", NORMAL);
        printn(gBestFit, "|");
        arr::print(gBest.data(), p.dims);
    }

    void
    showPopulation(const Population & pop,
                   const int popSize,
                   const char * const title)
    {
        print(BLUE, title, popSize, NORMAL);
        for (int i=0;i!=popSize;++i)
        {
            if (i % 2) printn(DARKER);

            auto & w = pop.weights[i];
            print(i, "|", w.pInertia, w.pMemory, w.pCooperation, w.pPerturbation, w.dThreshold, w.dVelocity);

            printn("    ");
            arr::print(&pop.particles(i,0), p.dims);

            printn("    ");
            arr::print(&pop.velocity(i,0), p.dims);
            if (i % 2) printn(NORMAL);
        }
    }

};

#endif // CDEEPSO_HPP
