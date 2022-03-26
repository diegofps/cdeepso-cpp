#ifndef CDEEPSO_PARAMS_HPP
#define CDEEPSO_PARAMS_HPP

#include <wup/wup.hpp>

using namespace wup;

typedef double Precision;
typedef wup::base_random<NaiveGenerator<Precision>, Precision> Random;

class CDEEPSOParams
{
public:

    enum MemStrategy {
        POS=1,
        MEM=2,
        POS_MEM=3
    };

    enum DEType {
        RAND=2,
        BEST=3
    };

private:

    class MemStrategyDecoder : public std::map<std::string, MemStrategy>
    {
    public:
        MemStrategyDecoder()
        {
            (*this)["POS"] = MemStrategy::POS;
            (*this)["MEM"] = MemStrategy::MEM;
            (*this)["POS_MEM"] = MemStrategy::POS_MEM;
        }
    };

    class DETypeDecoder : public std::map<std::string, DEType>
    {
    public:
        DETypeDecoder()
        {
            (*this)["RAND"] = DEType::RAND;
            (*this)["BEST"] = DEType::BEST;
        }
    };

public:

    MemStrategy memStrategy = MemStrategy::MEM;
    DEType deType = DEType::BEST;

    Precision mutationRate = 0.5;
    Precision communicationProbability = 0.1;
    Precision maxVelocity = 2.0;
    Precision xMin = -1.0;
    Precision xMax = 1.0;

//    int blockSize = 10;
    int dims = 50;
    int popSize = 50;
    int memGBestSize = 5;
    int maxFitEval = 100000;
    int maxGen = 50000;
    int maxGenWoChangeBest = 1000;
    int printConvergenceResults = 100;
    int maxRun = 50;
    int threads = 0;
    int ntupleDims = 8;

    std::string eval = "ras";

public:

    CDEEPSOParams()
    {

    }

    CDEEPSOParams(Params & params)
    {
        //parseParams(params.use("detector"));
        parseParams(params);
    }
    
    void
    parseParams(Params & p)
    {
        p.popEnum<MemStrategyDecoder>("memStrategy", memStrategy);
        p.popEnum<DETypeDecoder>("deType", deType);

        p.popDouble("mutationRate", mutationRate);
        p.popDouble("communicationProbability", communicationProbability);
        p.popDouble("maxVelocity", maxVelocity);
        p.popDouble("xMin", xMin);
        p.popDouble("xMax", xMax);

//        p.popInt("blockSize", blockSize);
        p.popInt("dims", dims);
        p.popInt("popSize", popSize);
        p.popInt("memGBestSize", memGBestSize);
        p.popInt("maxFitEval", maxFitEval);
        p.popInt("maxGen", maxGen);
        p.popInt("maxGenWoChangeBest", maxGenWoChangeBest);
        p.popInt("printConvergenceResults", printConvergenceResults);
        p.popInt("maxRun", maxRun);
        p.popInt("threads", threads);
        p.popInt("ntupleDims", ntupleDims);

        p.popString("eval", eval);
    }

    void
    display()
    {
        print(YELLOW, "\n--- CDEEPSO++ Parameters ---\n", WHITE);

        print("memStrategy =", memStrategy);
        print("deType =", deType);

        print("mutationRate =", mutationRate);
        print("communicationProbability =", communicationProbability);
        print("maxVelocity =", maxVelocity);

//        print("blockSize =", blockSize);
        print("dims =", dims);
        print("xMin =", xMin);
        print("xMax =", xMax);
        print("popSize =", popSize);
        print("memGBestSize =", memGBestSize);
        print("maxFitEval =", maxFitEval);
        print("maxGen =", maxGen);
        print("maxGenWoChangeBest =", maxGenWoChangeBest);
        print("printConvergenceResults =", printConvergenceResults);
        print("maxRun =", maxRun);
        print("threads =", threads);
        print("ntupleDims =", ntupleDims);

        print("eval =", eval);

        printn(NORMAL);
    }

};

#endif // CDEEPSO_PARAMS_HPP
