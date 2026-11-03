#include "AACodonMutSelOmegaCodonSubMatrix.hpp"
#include <tuple>

void AACodonMutSelOmegaCodonSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < Nstate; i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) * GetCodonFitness(i);
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < Nstate; i++) { mStationary[i] /= total; }
}

void AACodonMutSelOmegaCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (auto j : statespace->GetNeighbors(i)) {
        int pos = GetDifferingPosition(i, j);
        int a = GetCodonPosition(pos, i);
        int b = GetCodonPosition(pos, j);

        Q(i, j) = (*NucMatrix)(a, b);
        Q(i, j) *= FixationFactor(GetDeltaS(i, j));

        if (!Synonymous(i, j)) { Q(i, j) *= GetOmega(); }

        total += Q(i, j);

        assert(!std::isinf(Q(i, j)));
        assert(Q(i, j) >= 0);
    }

    Q(i, i) = -total;
    assert(total >= 0);
}

std::tuple<double, double> AACodonMutSelOmegaCodonSubMatrix::GetFlowDNDS() const {
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

            double flow = (*NucMatrix)(a, b) * FixationFactor(GetDeltaS(i, j));
            if (!Synonymous(i, j)) {
                dn += flow;
            } else {
                ds += flow;
            }
        }

        totdn += mStationary[i] * dn;
        totds += mStationary[i] * ds;
    }
    return std::make_tuple(totdn, totds);
}

std::tuple<double, double> AACodonMutSelOmegaCodonSubMatrix::GetRelativeFlowDNDS() const {
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
            double pfix = FixationFactor(GetDeltaS(i, j));

            if (!Synonymous(i, j)) {
                mutdn += nucrate;
                subdn += nucrate * pfix;
            } else {
                mutds += nucrate;
                subds += nucrate * pfix;
            }
        }

        totmutdn += mStationary[i] * mutdn;
        totmutds += mStationary[i] * mutds;
        totsubdn += mStationary[i] * subdn;
        totsubds += mStationary[i] * subds;
    }

    return std::make_tuple(totsubdn / totmutdn, totsubds / totmutds);
}

double AACodonMutSelOmegaCodonSubMatrix::GetPredictedDNDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetFlowDNDS();
    return dn / ds;
}

double AACodonMutSelOmegaCodonSubMatrix::GetPredictedDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetFlowDNDS();
    return ds;
}

double AACodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDNDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return dn / ds;
}

double AACodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDN() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return dn;
}

double AACodonMutSelOmegaCodonSubMatrix::GetPredictedRelativeDS() const {
    double dn = 0, ds = 0;
    std::tie(dn, ds) = GetRelativeFlowDNDS();
    return ds;
}
