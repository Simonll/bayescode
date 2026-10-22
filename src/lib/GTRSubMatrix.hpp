#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "SubMatrix.hpp"

/**
 * \brief A general time reversible substitution matrix
 *
 * Parameterized in terms of an array of relative rates (size
 * Nstate*(Nstate-1)/2) and an array of equilibrium frequencies (size Nstate)
 */

class GTRSubMatrix : public virtual SubMatrix {
  public:
    //! constructor parameterized by an array of relative rates (size
    //! Nstate*(Nstate-1)/2) and an array of equilibrium frequencies (size Nstate)
    GTRSubMatrix(int inNstate, const std::vector<double> &rr, const std::vector<double> &stat,
        bool innormalise = false);
    ~GTRSubMatrix() override = default;

    //! return number of relative rates
    int GetNRelativeRate() const { return Nrr; }
    //! converter returning relrate(i,j) == relrate(j,i), for any pair of states
    //! i!=j
    double RelativeRate(int i, int j) const { return mRelativeRate[rrindex(i, j)]; }

    //! make a copy of the entries of the equilibrium frequency vector; should be
    //! done each time this vector has been modified
    void CopyStationary(const std::vector<double> &instat);

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override {}
    double GetRate() const override;

    const std::vector<double> &mRelativeRate;
    int Nrr;

  private:
    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * Nstate - j - 1) * j / 2 + i - j - 1;
    }
};

//! GTRSubMatrix with the x3 codon factor in GetRate().
//!
//! The nucleotide GTR is normalised so that the neutral mutation flux is
//! 1 per codon (3 positions x 1/3 per nt). The x3 is only needed by the
//! codon mut-sel model (CodonMutSelMultipleOmega); the base GTRSubMatrix
//! keeps the plain 2*norm scale for the other models (SingleOmega, M2a, ...).
class GTRSubMatrixCodon : public GTRSubMatrix {
  public:
    GTRSubMatrixCodon(int inNstate, const std::vector<double> &rr,
        const std::vector<double> &stat, bool innormalise = false)
        : SubMatrix(inNstate, innormalise), GTRSubMatrix(inNstate, rr, stat, innormalise) {}

  protected:
    double GetRate() const override;
};

#endif
