#pragma once

#include <cassert>
#include <cmath>
#include <tuple>
#include "CodonMutSelOmegaCodonSubMatrix.hpp"  // FixationFactor, kFitnessFloor
#include "CodonSubMatrix.hpp"

/**
 * \brief A mutation-selection codon substitution process with amino-acid and codon fitnesses.
 *
 * This codon substitution process describes the evolution of a coding position,
 * under a constant fitness landscape combining a 20-vector of amino-acid fitnesses and a
 * 61-vector of codon fitnesses (codon preference). The model is parameterized by a nucleotide
 * matrix, specifying the mutation process, a vector of 20 scaled fitness parameters (summing to
 * 1) for the amino-acids, and a vector of 61 scaled fitness parameters for the codons.
 *
 * - stationary distribution: pi_nuc(codon) x fitness_AA(aa(codon)) x fitness_codon(codon);
 * - non-synonymous deltaS = dlogFitness(AA) + dlogCodonFitness;
 * - synonymous deltaS = dlogCodonFitness.
 *
 * The process also takes a real parameter, omega, which acts as a multiplier in front of all
 * non-synonymous substitutions. The standard mutation-selection process is obtained by setting
 * omega=1 (Rodrigue and Lartillot, 2017, MBE).
 */

class AACodonMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix,
                                         public virtual OmegaCodonSubMatrix {
  public:
    //! constructor, parameterized by a codon state space (genetic code), a
    //! nucleotide mutation matrix, a 61-vector of codon fitnesses, a 20-vector of amino-acid
    //! fitnesses, and a positive real parameter omega (=1 in the standard model).
    AACodonMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace,
        const SubMatrix *inNucMatrix, const std::vector<double> &incodon,
        const std::vector<double> &inaa, double inomega, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise),
          fitnesses(inaa.size(), 0.0),
          logfitnesses(inaa.size(), 0.0),
          codonfitnesses(incodon.size(), 0.0),
          logcodonfitnesses(incodon.size(), 0.0),
          aa(inaa),
          codon(incodon) {
        RecomputeFitness();
    }

    //! \brief access by copy to fitness of a given amino-acid
    //!
    //! Note: to avoid numerical errors, this function adds kFitnessFloor (1e-8).
    double GetFitness(int a) const {
        assert(std::abs((aa[a] + kFitnessFloor) - fitnesses[a]) < 1e-12);
        return fitnesses[a];
    }

    double GetLogFitness(int a) const {
        assert(std::abs(log(GetFitness(a)) - logfitnesses[a]) < 1e-12);
        return logfitnesses[a];
    }

    //! access to the fitness of a given codon (with kFitnessFloor added)
    double GetCodonFitness(int c) const {
        assert(std::abs((codon[c] + kFitnessFloor) - codonfitnesses[c]) < 1e-12);
        return codonfitnesses[c];
    }

    double GetLogCodonFitness(int c) const {
        assert(std::abs(log(GetCodonFitness(c)) - logcodonfitnesses[c]) < 1e-12);
        return logcodonfitnesses[c];
    }

    std::tuple<double, double> GetFlowDNDS() const;
    std::tuple<double, double> GetRelativeFlowDNDS() const;

    double GetPredictedDNDS() const;
    double GetPredictedDS() const;

    double GetPredictedRelativeDNDS() const;
    double GetPredictedRelativeDN() const;
    double GetPredictedRelativeDS() const;

    void CorruptMatrixNoFitnessRecomput() { SubMatrix::CorruptMatrix(); }

    void CorruptMatrix() override {
        RecomputeFitness();
        SubMatrix::CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    //! deltaS (log space) of the substitution i -> j: dlogFitness(AA) + dlogCodonFitness
    //! for non-synonymous changes, dlogCodonFitness only for synonymous ones.
    double GetDeltaS(int i, int j) const {
        double deltaS = GetLogCodonFitness(j) - GetLogCodonFitness(i);
        if (!Synonymous(i, j)) {
            deltaS += GetLogFitness(GetCodonStateSpace()->Translation(j)) -
                      GetLogFitness(GetCodonStateSpace()->Translation(i));
        }
        return deltaS;
    }

    //! recompute fitnesses and log-fitnesses from the aa and codon vectors (with floor)
    void RecomputeFitness() {
        for (size_t a{0}; a < aa.size(); a++) {
            fitnesses[a] = aa[a] + kFitnessFloor;
            logfitnesses[a] = log(fitnesses[a]);
        }
        for (size_t c{0}; c < codon.size(); c++) {
            codonfitnesses[c] = codon[c] + kFitnessFloor;
            logcodonfitnesses[c] = log(codonfitnesses[c]);
        }
    }

    // fitness precomputation
    std::vector<double> fitnesses;
    std::vector<double> logfitnesses;
    std::vector<double> codonfitnesses;
    std::vector<double> logcodonfitnesses;

    // data members
    const std::vector<double> &aa;
    const std::vector<double> &codon;
};
