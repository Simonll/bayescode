#include "GTRSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     GTRSubMatrix
// ---------------------------------------------------------------------------

GTRSubMatrix::GTRSubMatrix(
    int inNstate, const std::vector<double> &rr, const std::vector<double> &stat, bool innormalise)
    : SubMatrix(inNstate, innormalise), mRelativeRate(rr) {
    Nrr = Nstate * (Nstate - 1) / 2;
    CopyStationary(stat);
}

void GTRSubMatrix::CopyStationary(const std::vector<double> &instat) {
    for (int k = 0; k < Nstate; k++) { mStationary[k] = instat[k]; }
}

// ---------------------------------------------------------------------------
//     ComputeArray
// ---------------------------------------------------------------------------

void GTRSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < Nstate; j++) {
        if (i != j) {
            Q(i, j) = RelativeRate(i, j) * mStationary[j];
            total += Q(i, j);
        }
    }

    Q(i, i) = -total;
}

// ---------------------------------------------------------------------------
//     ComputeRate()
// ---------------------------------------------------------------------------

double GTRSubMatrix::GetRate() const {
    if (!ArrayUpdated()) {
        UpdateStationary();
        for (int k = 0; k < Nstate; k++) { ComputeArray(k); }
    }
    for (int k = 0; k < Nstate; k++) { flagarray[k] = true; }
    double norm = 0;
    for (int i = 0; i < Nstate - 1; i++) {
        for (int j = i + 1; j < Nstate; j++) { norm += mStationary[i] * Q(i, j); }
    }
    return 2 * norm;
}
// ---------------------------------------------------------------------------
//     GTRSubMatrixCodon: GetRate() with the x3 codon factor
// ---------------------------------------------------------------------------
// The nucleotide GTR is normalised so that the neutral mutation flux is
// 1 per codon (3 positions x 1/3 per nt). The x3 is only needed by the
// codon mut-sel model (CodonMutSelMultipleOmega); the base GTRSubMatrix
// keeps the plain 2*norm scale so that the other models (SingleOmega,
// M2a, ...) are unaffected.
double GTRSubMatrixCodon::GetRate() const {
    if (!ArrayUpdated()) {
        UpdateStationary();
        for (int k = 0; k < Nstate; k++) { ComputeArray(k); }
    }
    for (int k = 0; k < Nstate; k++) { flagarray[k] = true; }
    double norm = 0;
    for (int i = 0; i < Nstate - 1; i++) {
        for (int j = i + 1; j < Nstate; j++) { norm += mStationary[i] * Q(i, j); }
    }
    return 2 * norm * 3;
}
