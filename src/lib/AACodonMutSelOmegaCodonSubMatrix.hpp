#pragma once

#include <cassert>
#include "CodonSubMatrix.hpp"

/**
 * \brief A mutation-selection codon substitution process.
 *
 * This codon substitution process describes the evolution of a coding position,
 * under a constant fitness landscape over the 20 amino-acids. The model is
 * parameterized by a nucleotide matrix, specifying the mutation process, and a
 * vector of 20 scaled fitness parameters (summing to 1), for the 20
 * amino-acids. The process also takes a real parameter, omega, which acts as a
 * multiplier in front of all non-synonymous substitutions. The standard
 * mutation-selection process is obtained by setting omega=1. Letting omega be
 * different from 1 was explored in Rodrigue and Lartillot, 2017, MBE (detecting
 * deviations from expected non-syn rate under the standard mut-sel model).
 */

class AACodonMutSelOmegaCodonSubMatrix : public virtual NucCodonSubMatrix,
                                         public virtual OmegaCodonSubMatrix {
  public:
    //! constructor, parameterized by a codon state space (genetic code), a
    //! nucleotide mutation matrix, a 20-vector of amino-acid fitnesss, and a
    //! positive real parameter omega (=1 in the standard model).
    AACodonMutSelOmegaCodonSubMatrix(const CodonStateSpace *instatespace,
        const SubMatrix *inNucMatrix, const std::vector<double> &incodon,
        const std::vector<double> &inaa, double inomega, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          OmegaCodonSubMatrix(instatespace, inomega, innormalise),
          codonfitnesses(incodon.size(), 0.0),
          logcodonfitnesses(incodon.size(), 1.0 / incodon.size()),
          fitnesses(inaa.size(), 0.0),
          logfitnesses(inaa.size(), 1.0 / inaa.size()),
          aa(inaa),
          codon(incodon) {}

    //! \brief access by copy to fitness of a given amino-acid
    //!
    //! Note: to avoid numerical errors, this function adds 1e-8.
    double GetFitness(int a) const {
        assert(std::abs((exp(log(aa[a])) + 1e-8) - fitnesses[a]) < 1e-6);
        return fitnesses[a];
    }

    double GetLogFitness(int a) const {
        assert(std::abs(log(GetFitness(a)) - logfitnesses[a]) < 1e-6);
        return logfitnesses[a];
    }

    double GetCodonFitness(int a) const {
        assert(std::abs((exp(log(codon[a])) + 1e-8) - codonfitnesses[a]) < 1e-6);
        return codonfitnesses[a];
    }

    double GetLogCodonFitness(int a) const {
        assert(std::abs(log(GetCodonFitness(a)) - logcodonfitnesses[a]) < 1e-6);
        return logcodonfitnesses[a];
    }

    std::tuple<double, double> GetFlowDNDS() const;
    std::tuple<double, double> GetRelativeFlowDNDS() const;

    double GetPredictedDNDS() const;

    double GetPredictedDS() const;

    double GetPredictedRelativeDNDS() const;

    double GetPredictedRelativeDS() const;

    void CorruptMatrixNoFitnessRecomput() { SubMatrix::CorruptMatrix(); }

    void CorruptMatrix() override {
        for (size_t a{0}; a < aa.size(); a++) {
            fitnesses[a] = exp(log(aa[a])) + 1e-8;
            logfitnesses[a] = log(fitnesses[a]);
        }
        for (size_t c{0}; c < codon.size(); c++) {
            codonfitnesses[c] = exp(log(codon[c])) + 1e-8;
            logcodonfitnesses[c] = log(codonfitnesses[c]);
        }
        SubMatrix::CorruptMatrix();
    }

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    // fitness precomputation
    std::vector<double> fitnesses;
    std::vector<double> logfitnesses;
    std::vector<double> codonfitnesses;
    std::vector<double> logcodonfitnesses;

    // data members
    const std::vector<double> &aa;
    const std::vector<double> &codon;
};