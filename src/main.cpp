#include "cdeepso.hpp"
#include "cdeepso_params.hpp"
#include "functions.hpp"

#include <iostream>
#include <wup/wup.hpp>
#include <sstream>

WUP_STATICS;

using namespace std;
using namespace wup;

void
custom(Population & pop,
       Refreshes & refresh,
       Fitness & fitness)
{
    for (uint i=0;i!=pop.particles.numRows();++i)
        if (refresh[i])
            fitness[i] = rosenbrock(&pop.particles(i,0), pop.particles.numCols());
}

int
main(const int argc, const char * argv[])
{
    srand(time(NULL));
    printn(std::scientific);

    Params p(argc, argv);
    CDEEPSOParams cp(p);
    vector<Precision> allFits(cp.maxRun);
    vector<long double> ellapsed(cp.maxRun);

    void (*eval)(Particles & particles,
                 Refreshes & refresh,
                 Fitness & fitness) = nullptr;

    if (cp.eval == "ras") eval = rastrigin;
    else if (cp.eval == "ros") eval = rosenbrock;
    else if (cp.eval == "gri") eval = griewank;
    else error("Invalid eval function:", cp.eval);

    cp.display();

    Clock cc;

    print(YELLOW, "\n--- CDEEPSO++ Main Loop ---\n", NORMAL);

    if (cp.threads == 1)
    {
        Clock c;
        for (int r=0;r!=cp.maxRun;++r)
        {
            c.start();
            CDEEPSO m(cp);

            m.optimize(eval);

            ellapsed[r] = c.lap_milli();
            allFits[r] = m.gBestFit;
//            printn(cat(WHITE, r, NORMAL, " : Best fitness = ", GREEN, m.gBestFit, "\n", NORMAL));
        }
    }

    else
    {
        wup::parallel(cp.threads, cp.maxRun, [&](const int tid, const int jid) {
            UNUSED(tid);

            Clock c;
            CDEEPSO m(cp);

            m.optimize(eval);

            ellapsed[jid] = c.stop().ellapsed_milli();
            allFits[jid] = m.gBestFit;
//            printn(cat(WHITE, jid, NORMAL, " : Best fitness = ", GREEN, m.gBestFit, "\n", NORMAL));
        });
    }

    long double totalTime = cc.stop().ellapsed_milli();
    long double minimum, maximum, mean, std;

    print(YELLOW, "\n--- CDEEPSO++ Results ---\n", WHITE);

    arr::stats(allFits, minimum, maximum, mean, std);
    print("Fitness:");
    print("  Minimum:", minimum);
    print("  Maximum:", maximum);
    print("  Mean:", mean);
    print("  Std:", std);

    print("Total execution time:", totalTime, "ms");

    arr::stats(ellapsed, minimum, maximum, mean, std);
    print("Time to run:");
    print("  Minimum:", minimum, "ms");
    print("  Maximum:", maximum, "ms");
    print("  Mean:", mean, "ms");
    print("  Std:", std, "ms");
    printn(NORMAL);

    return 0;
}
