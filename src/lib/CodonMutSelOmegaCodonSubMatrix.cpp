#include "CodonMutSelOmegaCodonSubMatrix.hpp"
#include <tuple>

void CodonMutSelOmegaCodonSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < Nstate; i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) * GetFitness(i);
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    // double min = 1;
    for (int i = 0; i < Nstate; i++) { mStationary[i] /= total; }
}

void CodonMutSelOmegaCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (auto j : statespace->GetNeighbors(i)) {
        int pos = GetDifferingPosition(i, j);
        int a = GetCodonPosition(pos, i);
        int b = GetCodonPosition(pos, j);

        Q(i, j) = (*NucMatrix)(a, b);

        double deltaS = GetLogFitness(j) - GetLogFitness(i);
        Q(i, j) *= FixationFactor(deltaS);

        if (!Synonymous(i, j)) { Q(i, j) *= GetOmega(); }

        total += Q(i, j);

        assert(!std::isinf(Q(i, j)));
        assert(Q(i, j) >= 0);
    }

    Q(i, i) = -total;
    assert(total >= 0);
}

std::tuple<double, double> CodonMutSelOmegaCodonSubMatrix::GetFlowDNDS() const {
    GetStationary();  // recomputed only if the matrix was corrupted since last use
    double totdn = 0;
    double totds = 0;
    for (int i = 0; i < Nstate; i++) {
        double ds = 0;
        double dn = 0;
        for (auto j : statespace->GetNeighbors(i)) {
            int pos = GetDifferingPosition(i, j);
            int a = GetCodonPosition(pos, i);
            int b = GetCodonPosition(pos, j);

            double nucrate = (*NucMatrix)(a, b);
            double deltaS = GetLogFitness(j) - GetLogFitness(i);
            double pfix = FixationFactor(deltaS);
            if (!Synonymous(i, j)) {
            
                dn += nucrate * pfix;

            } else {
               
                ds += nucrate * pfix;
            }
        }

        totdn += mStationary[i] * dn;
        totds += mStationary[i] * ds;
    }
    return std::make_tuple(totdn, totds);
}

double CodonMutSelOmegaCodonSubMatrix::GetPredictedDNDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetFlowDNDS();
    return dn / ds;
}

double CodonMutSelOmegaCodonSubMatrix::GetPredictedDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetFlowDNDS();
    return ds;
}

double CodonMutSelOmegaCodonSubMatrix::GetPredictedDN() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetFlowDNDS();
    return dn;
}

std::tuple<double, double> CodonMutSelOmegaCodonSubMatrix::GetRelativeFlowDNDS() const {
    GetStationary();  // recomputed only if the matrix was corrupted since last use

    double totmutdn = 0;
    double totmutds = 0;

    double totsubdn = 0;
    double totsubds = 0;

    for (int i = 0; i < Nstate; i++) {
        double mutds = 0;
        double mutdn = 0;
        double subds = 0;
        double subdn = 0;
        for (auto j : statespace->GetNeighbors(i)) {
            int pos = GetDifferingPosition(i, j);
            int a = GetCodonPosition(pos, i);
            int b = GetCodonPosition(pos, j);

            double nucrate = (*NucMatrix)(a, b);
            double deltaS = GetLogFitness(j) - GetLogFitness(i);
            double pfix = FixationFactor(deltaS);
            if (!Synonymous(i, j)) {

                mutdn += nucrate;
                subdn += nucrate * pfix;

            } else {
               
                mutds += nucrate;
                subds += nucrate * pfix;
            }
        }

        totsubdn += mStationary[i] * subdn;
        totmutdn += mStationary[i] * mutdn;       
        totsubds += mStationary[i] * subds;
        totmutds += mStationary[i] * mutds;
    }

    // guard against division by zero (degenerate state space)
    double rdn = (totmutdn > 0) ? totsubdn / totmutdn : 0.0;
    double rds = (totmutds > 0) ? totsubds / totmutds : 0.0;
    return std::make_tuple(rdn, rds);
}


double CodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDNDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return dn / ds;
}

double CodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return ds;
}

double CodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDN() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return dn;
}