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
        if ((fabs(deltaS)) < 1e-30) {
            Q(i, j) *= 1 + deltaS / 2;
        } else if (deltaS > 50) {
            Q(i, j) *= deltaS;
        } else if (deltaS < -50) {
            Q(i, j) = 0;
        } else {
            Q(i, j) *= deltaS / (1.0 - exp(-deltaS));
        }

        if (!Synonymous(i, j)) { Q(i, j) *= GetOmega(); }

        total += Q(i, j);

        assert(!std::isinf(Q(i, j)));
        assert(Q(i, j) >= 0);
    }

    Q(i, i) = -total;
    assert(total >= 0);
}

std::tuple<double, double> CodonMutSelOmegaCodonSubMatrix::GetFlowDNDS() const {
    UpdateStationary();
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

            if (!Synonymous(i, j)) {
                double deltaS = GetLogFitness(j) - GetLogFitness(i);
                double pfix;
                if ((fabs(deltaS)) < 1e-30) {
                    pfix = 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    pfix = deltaS;
                } else if (deltaS < -50) {
                    pfix = 0;
                } else {
                    pfix = deltaS / (1.0 - exp(-deltaS));
                }

                dn += nucrate * pfix;

            } else {
                double deltaS = GetLogFitness(j) - GetLogFitness(i);
                double pfix;
                if ((fabs(deltaS)) < 1e-30) {
                    pfix = 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    pfix = deltaS;
                } else if (deltaS < -50) {
                    pfix = 0;
                } else {
                    pfix = deltaS / (1.0 - exp(-deltaS));
                }

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


std::tuple<double, double> CodonMutSelOmegaCodonSubMatrix::GetRelativeFlowDNDS() const {
    UpdateStationary();


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

            if (!Synonymous(i, j)) {
                double deltaS = GetLogFitness(j) - GetLogFitness(i);
                double pfix;
                if ((fabs(deltaS)) < 1e-30) {
                    pfix = 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    pfix = deltaS;
                } else if (deltaS < -50) {
                    pfix = 0;
                } else {
                    pfix = deltaS / (1.0 - exp(-deltaS));
                }

                mutdn += nucrate;
                subdn += nucrate * pfix;

            } else {
                double deltaS = GetLogFitness(j) - GetLogFitness(i);
                double pfix;
                if ((fabs(deltaS)) < 1e-30) {
                    pfix = 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    pfix = deltaS;
                } else if (deltaS < -50) {
                    pfix = 0;
                } else {
                    pfix = deltaS / (1.0 - exp(-deltaS));
                }
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