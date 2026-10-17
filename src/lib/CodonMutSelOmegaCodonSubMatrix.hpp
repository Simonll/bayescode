#pragma once

#include <cassert>
#include <cmath>
#include "CodonSubMatrix.hpp"

/**
 * \brief A mutation-selection codon substitution process.
 *
 * This codon substitution process describes the evolution of a coding position,
 * under a constant fitness landscape over the codon state space, including amino acids and codons
 * preferences. The model is parameterized by a nucleotide matrix, specifying the mutation process,
 * and a vector of 61 scaled fitness parameters (summing to 1), for the 61 codons. The process
 * also takes a real parameter, omega, which acts as a multiplier in front of all non-synonymous
 * substitutions. The standard mutation-selection process is obtained by setting omega=1. Letting
 * omega be different from 1 was explored in Rodrigue and Lartillot, 2017, MBE (detecting deviations
 * from expected non-syn rate under the standard mut-sel model).
 */

//! fixation probability factor for a fitness difference deltaS (in log space).
//! Shared by ComputeArray, GetFlowDNDS, GetRelativeFlowDNDS and the reader, so a
//! change to the biology is applied in one place only.
inline double FixationFactor(double deltaS) {
    if (std::fabs(deltaS) < 1e-8) { return 1 + deltaS / 2; }
    if (deltaS > 50) { return deltaS; }
    if (deltaS < -50) { return 0; }
    return deltaS / -std::expm1(-deltaS);
}

//! floor added to codon fitnesses to avoid log(0) on (near-)zero fitnesses.
//! Shared by the matrix and the reader so the two agree on the fitness values.
constexpr double kFitnessFloor = 1e-8;

class CodonMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix,
                                       public virtual OmegaCodonSubMatrix {
  public:
    //! constructor, parameterized by a codon state space (genetic code), a
    //! nucleotide mutation matrix, a 61-vector of codon fitnesses, and a
    //! positive real parameter omega (=1 in the standard model).
    CodonMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace,
        const SubMatrix *inNucMatrix, const std::vector<double> &incodon, double inomega,
        bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise),
          fitnesses(incodon.size(), 0.0),
          logfitnesses(incodon.size(), 0.0),
          codon(incodon) {
        RecomputeFitness();
    }

    //! \brief access by copy to fitness of a given codon
    //!
    //! Note: to avoid numerical errors, this function adds 1e-8.
    double GetFitness(int c) const {
        assert(std::abs((codon[c] + kFitnessFloor) - fitnesses[c]) < 1e-6);
        return fitnesses[c];
    }

    double GetLogFitness(int c) const {
        assert(std::abs(log(GetFitness(c)) - logfitnesses[c]) < 1e-6);
        return logfitnesses[c];
    }

    std::tuple<double, double> GetFlowDNDS() const;
    double GetPredictedDNDS() const;
    double GetPredictedDS() const;
    double GetPredictedDN() const;

    std::tuple<double, double> GetRelativeFlowDNDS() const;

    double GetPredictedRelativeDNDS() const;
    double GetPredictedRelativeDS() const;
    double GetPredictedRelativeDN() const;

    void CorruptMatrixNoFitnessRecomput() { SubMatrix::CorruptMatrix(); }

    void CorruptMatrix() override {
        RecomputeFitness();
        SubMatrix::CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    //! \brief normalisation scale for the codon matrix.
    //!
    //! Returns the nucleotide mutation scale (independent of omega), matching
    //! pbmpi's CodonMutSelProfileSubMatrix::GetRate. Normalising the codon
    //! matrix by this scale (rather than by its own flux) keeps omega as a pure
    //! syn/non-syn ratio parameter: it does not change the global rate of the
    //! matrix, and the branch lengths stay comparable across datasets.
    double GetRate() const override;

    //! recompute fitnesses and log-fitnesses from the codon vector (with floor)
    void RecomputeFitness() {
        for (size_t c{0}; c < codon.size(); c++) {
            fitnesses[c] = codon[c] + kFitnessFloor;
            logfitnesses[c] = log(fitnesses[c]);
        }
    }
    // fitness precomputation
    std::vector<double> fitnesses;
    std::vector<double> logfitnesses;

    // data members
    const std::vector<double> &codon;
};