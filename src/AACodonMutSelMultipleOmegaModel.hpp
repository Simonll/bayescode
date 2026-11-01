#pragma once

#include <algorithm>
#include <limits>
#include "AcceptanceTracker.hpp"
#include "AACodonMutSelCodonMatrixBidimArray.hpp"
#include "Chrono.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "StickBreakingProcess.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"
#include "components/model_decl_utils.hpp"
#include "lib/Dirichlet.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- double Dirichlet process version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and
 * amino-acid a
 * - an array of site-specific omega multiplier, capturing deviations of the non-syn rate from the
 * model; this parameter is fixed to 1 by default.
 * - an array of codon fitness prfole F_c, for codon c
 *
 * Site-specific amino-acid fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration
 * parameter kappa, and truncated at Ncat.
 *
 * The base distribution of this Dirichlet process, G_0, is itself a (truncated
 * stick breaking) mixture of Dirichlet distributions, parameterized by
 * basekappa, and truncated at baseNcat. Each component of this mixture
 * Dirichlet is parameterized by a center (a 20-freqeuncy vector) and a
 * concentration.
 *
 * Concerning the first-level stick-breaking process, by default, Ncat == 100
 * (can be changed using the  -ncat option). As for baseNcat, it is equal to 1,
 * in which case the base distribution G_0 is just a simple Dirichlet (as in
 * Rodrigue et al, 2010). This simple model is probably the one that should be
 * used by default for single-gene analyses. The more complex model with
 * baseNcat > 1 is meant to be used in a multi-gene context (although, even in
 * that case, mixing is still challenging, and not sure whether this model
 * brings important improvement in the end). Thus, baseNcat = 1 appears to be
 * the most reasonable model settings for now.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - codonfinesses uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 20.
 *
 */


// Input-file parsing: same strict behaviour as in CodonMutSelMultipleOmegaModel.hpp (fatal error
// with a message on malformed input, rejection of empty files and non-positive values). The
// helpers live in their own namespace so that this header can coexist in a translation unit with
// the other model headers, which define their own (differently-behaving) copies at global scope.
namespace aacodonmutsel {

inline double distance(std::vector<double> const &v1, std::vector<double> const &v2) {
    double tot = 0;
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); i++) { tot += std::abs(v1[i] - v2[i]); }
    return tot;
}

//! Fatal error while parsing an input file
[[noreturn]] inline void parsing_error(std::string const &file_name, std::string const &what) {
    std::cerr << "error in file " << file_name << ": " << what << std::endl;
    exit(1);
}

//! Guess the separator (space, comma or tab) from the header line
inline char detect_separator(std::string const &header) {
    char sep{' '};
    long best = 0;
    for (char sep_test : std::vector<char>({' ', ',', '\t'})) {
        long n = std::count(header.begin(), header.end(), sep_test);
        if (n > best) {
            sep = sep_test;
            best = n;
        }
    }
    return sep;
}

//! Split a line on sep, keeping empty inner fields; a single trailing separator is ignored
inline std::vector<std::string> split_fields(std::string line, char sep) {
    if (!line.empty() and line.back() == '\r') { line.pop_back(); }
    if (!line.empty() and line.back() == sep) { line.pop_back(); }
    std::vector<std::string> fields;
    std::string word;
    std::istringstream line_stream(line);
    while (getline(line_stream, word, sep)) { fields.push_back(word); }
    if (!line.empty() and line.back() == sep) { fields.emplace_back(); }
    return fields;
}

//! Strict conversion: the whole field must be a finite number
inline double parse_number(std::string const &file_name, std::string const &word) {
    size_t pos = 0;
    double value = 0;
    try {
        value = std::stod(word, &pos);
    } catch (std::exception const &) { pos = 0; }
    while (pos > 0 and pos < word.size() and std::isspace(static_cast<unsigned char>(word[pos]))) {
        pos++;
    }
    if (pos == 0 or pos != word.size() or !std::isfinite(value)) {
        parsing_error(file_name, "'" + word + "' is not a valid number");
    }
    return value;
}

//! Read delta-omega values, one per line; each must be a finite number > 0
inline std::vector<double> open_delta_omega_array(std::string const &file_name) {
    std::vector<double> delta_omega_array{};
    std::ifstream input_stream(file_name);
    if (!input_stream) {
        std::cerr << "Delta-omega array file " << file_name << " doesn't exist" << std::endl;
        exit(1);
    }

    std::string line;
    size_t line_number = 0;
    while (getline(input_stream, line)) {
        line_number++;
        if (line.empty() or line == "\r") { continue; }
        double value = parse_number(file_name, line);
        // omega = 0 gives log(0) and a division by omega in GetOmega()
        if (!(value > 0.0)) {
            parsing_error(file_name, "line " + std::to_string(line_number) + ": omega value '" +
                                         line + "' must be strictly positive");
        }
        delta_omega_array.push_back(value);
    }
    if (delta_omega_array.empty()) {
        parsing_error(file_name, "no delta-omega values found; omegaNcat would be 0");
    }
    return delta_omega_array;
}

/**
 * Read the last `size` fields of each non-empty line after the header. Leading columns
 * (row names, ...) are ignored. All rows must have the same number of fields, at least `size`.
 */
inline std::vector<std::vector<double>> read_numeric_rows(
    std::string const &file_name, int size) {
    std::ifstream input_stream(file_name);
    if (!input_stream) {
        std::cerr << "file " << file_name << " doesn't exist" << std::endl;
        exit(1);
    }
    std::string line;
    if (!getline(input_stream, line)) { parsing_error(file_name, "file is empty"); }
    char sep = detect_separator(line);

    std::vector<std::vector<double>> rows;
    size_t nfields = 0;
    size_t line_number = 1;
    while (getline(input_stream, line)) {
        line_number++;
        if (line.empty() or line == "\r") { continue; }
        auto fields = split_fields(line, sep);
        if (rows.empty()) {
            nfields = fields.size();
            if (nfields < static_cast<size_t>(size)) {
                parsing_error(file_name, "line " + std::to_string(line_number) + " has " +
                                             std::to_string(nfields) +
                                             " fields, expected at least " + std::to_string(size));
            }
        } else if (fields.size() != nfields) {
            parsing_error(file_name, "line " + std::to_string(line_number) + " has " +
                                         std::to_string(fields.size()) + " fields, expected " +
                                         std::to_string(nfields));
        }
        std::vector<double> row(size);
        for (int i = 0; i < size; i++) {
            row[i] = parse_number(file_name, fields[nfields - size + i]);
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

//! Read amino-acid fitness profiles (`Nstate` columns, one profile per line after the header).
//! The first line of the file is always treated as a header and skipped; each profile must sum
//! to 1 and be strictly positive. Returns the distinct profiles and the allocation of each line.
inline std::tuple<std::vector<std::vector<double>>, std::vector<size_t>> open_preferences(
    std::string const &file_name, int Nstate) {
    std::vector<std::vector<double>> fitness_profiles{};
    std::vector<size_t> alloc;

    for (auto &fitness_profil : read_numeric_rows(file_name, Nstate)) {
        if (std::abs(std::accumulate(fitness_profil.begin(), fitness_profil.end(), 0.0) - 1.0) >
            1e-5) {
            parsing_error(file_name,
                "fitness profile " + std::to_string(alloc.size() + 1) + " doesn't sum to 1");
        }
        for (double fitness : fitness_profil) {
            // log(fitness + 1e-8) is NaN for a negative fitness, and a null one is degenerate
            if (!(fitness > 0.0)) {
                parsing_error(file_name, "fitness profile " + std::to_string(alloc.size() + 1) +
                                             " contains a non-positive value (" +
                                             std::to_string(fitness) + ")");
            }
        }
        bool push = true;
        for (size_t i = 0; i < fitness_profiles.size(); i++) {
            if (distance(fitness_profiles[i], fitness_profil) < 1e-5) {
                push = false;
                alloc.push_back(i);
                break;
            }
        }
        if (push) {
            alloc.push_back(fitness_profiles.size());
            fitness_profiles.push_back(std::move(fitness_profil));
        }
        assert(alloc[alloc.size() - 1] < fitness_profiles.size());
    }
    if (alloc.empty()) { parsing_error(file_name, "no preference profile found"); }
    return std::make_tuple(fitness_profiles, alloc);
}

/**
 * Read a single vector of `size` strictly positive values (header line, then one data line;
 * leading columns are ignored). The sum must be 1 up to 1e-5; the vector is then renormalized
 * to sum exactly to 1.
 */
inline std::vector<double> open_simplex_vector(
    std::string const &file_name, int size, std::string const &what) {
    auto rows = read_numeric_rows(file_name, size);
    if (rows.size() != 1) {
        parsing_error(file_name, what + ": expected exactly one data line after the header, found " +
                                     std::to_string(rows.size()));
    }
    std::vector<double> v = std::move(rows[0]);
    for (double x : v) {
        if (!(x > 0)) { parsing_error(file_name, what + ": values must be strictly positive"); }
    }
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    if (std::abs(sum - 1.0) > 1e-5) {
        parsing_error(file_name, what + " doesn't sum to 1 (sum = " + std::to_string(sum) + ")");
    }
    for (double &x : v) { x /= sum; }
    return v;
}

//! Read the vector of Nstate (61) codon fitnesses
inline std::vector<double> open_codonpreferences(std::string const &file_name, int Nstate) {
    return open_simplex_vector(file_name, Nstate, "codon fitness");
}

inline std::vector<double> open_nucstat(std::string const &file_name, int size) {
    return open_simplex_vector(file_name, size, "nucstat");
}

inline std::vector<double> open_nucrelrate(std::string const &file_name, int size) {
    return open_simplex_vector(file_name, size, "nucrelrate");
}

}  // namespace aacodonmutsel

class AACodonMutSelMultipleOmega : public ChainComponent {
    std::string datafile, treefile, aaprofilesfile{"Null"}, codonfitnessfile{"Null"},
        nucstatfile{"Null"}, nucrelratefile{"Null"};
    bool clamp_profiles{false}, clamp_profiles_allocation{false}, clamp_nucstat{false},
        clamp_nucrelrate{false}, clamp_codonfitness{false};

    std::unique_ptr<Tracer> tracer;
    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;
    int Nstate;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    SimpleBranchArray<double> *blhypermeanarray;
    GammaWhiteNoise *branchlength;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // nucleotide rates hyperparameters
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrixCodon *nucmatrix;

    // codonfitness hyperparameters
    std::vector<double> codonfitnesshypercenter;
    double codonfitnesshyperinvconc;

    // codonfitness
    std::vector<double> codonfitness;

    // SENCA-hybrid codon layer: `codonfitness` holds theta_c, the relative preference of codon c
    // among the synonymous codons of its amino acid (sum over syn(a) of theta_c = 1 for every a).
    // The substitution matrices use codonweight_c = n_a * theta_c (n_a = degeneracy of a), whose
    // mean is 1 within each AA: flat theta (theta_c = 1/n_a) <=> codonweight = 1 <=> plain AA
    // mutation-selection model. codonweight is derived from codonfitness (UpdateCodonWeights).
    std::vector<double> codonweight;
    // codon indices of each synonymous group (one group per AA present in the code)
    std::vector<std::vector<int>> synonymousgroups;

    // MCMC iteration counter (for periodic acceptance reporting)
    int accept_iter{0};

    bool clamp_delta_omega_array{false};
    std::string delta_omega_array_file = "Null";
    int omegaNcat;
    MultinomialAllocationVector *omega_alloc;

    double omega_weight_kappa;
    Dirichlet *omega_weight;
    double omega_shift;
    double delta_omegahypermean;
    double delta_omegahyperinvshape;
    // delta_omega_array of mean delta_omegahypermean and inverse shape parameter
    // delta_omegahyperinvshape
    IIDGamma *delta_omega_array;

    OmegaPathSuffStatArray *siteomegapathsuffstatarray;
    OmegaPathSuffStatArray *componentomegapathsuffstatarray;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions
    int baseNcat;
    int basemin;
    double basekappa;
    StickBreakingProcess *baseweight;

    OccupancySuffStat *baseprofile_occupancy;
    std::vector<double> basecenterhypercenter;
    double basecenterhyperinvconc;

    IIDDirichlet *basecenterarray;
    double baseconchypermean;
    double baseconchyperinvshape;

    IIDGamma *baseconcentrationarray;
    MultinomialAllocationVector *componentalloc;
    MixtureSelector<std::vector<double>> *componentcenterarray;
    MixtureSelector<double> *componentconcentrationarray;

    // aa fitness arrays across sites are a SBDP process of base G0 defined above
    int Ncat;
    double kappa;
    StickBreakingProcess *weight;

    OccupancySuffStat *profile_occupancy;
    OccupancySuffStat *omega_occupancy;

    MultiDirichlet *componentaafitnessarray;
    DirichletSuffStatArray *basesuffstatarray;

    MultinomialAllocationVector *profile_alloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
    AACodonMutSelCodonMatrixBidimArray *componentcodonmatrixbidimarray;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;
    DoubleMixtureSelector<AACodonMutSelOmegaCodonSubMatrix> *sitecodonsubmatrixarray;

    MixtureSelector<std::vector<double>> *siteaafitnessarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: free wo shrinkage
    int blmode;
    // currently, free wo shrinkages
    int nucmode;
    // currently, free wo shrinkage
    int basemode;

    // currently: free wo shrinkage or fixed (1 or 3)flatcodonfitness
    int omegamode;

    bool flataafitness;
    bool flatcodonfitness;
    bool flatnucstat;
    bool flatnucrelrate{false};  // no constructor argument: the nuc relative rates are free

    Chrono aachrono;
    Chrono basechrono;
    Chrono totchrono;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - inomegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    //! - Ncat: truncation of the first-level stick-breaking process (by default:
    //! 100)
    //! - baseNcat: truncation of the second-level stick-breaking process (by
    //! default: 1)
    AACodonMutSelMultipleOmega(std::string indatafile, std::string intreefile,
        std::string inaaprofilesfile, std::string incodonprofilesfile, std::string innucstat,
        std::string innucrelrate, int inomegamode, int inNcat, int inbaseNcat, int inomegaNcat,
        double inomegashift, bool inflatfitness, bool inflatcodonfitness,
        std::string indelta_omega_array_file, bool inflatnucstat, bool inflatnucrelrate = false)
        : datafile(std::move(indatafile)),
          treefile(std::move(intreefile)),
          aaprofilesfile(std::move(inaaprofilesfile)),
          codonfitnessfile(std::move(incodonprofilesfile)),
          nucstatfile(std::move(innucstat)),
          nucrelratefile(std::move(innucrelrate)),
          delta_omega_array_file(std::move(indelta_omega_array_file)),
          omegaNcat(inomegaNcat),
          omega_shift(inomegashift),
          baseNcat(inbaseNcat),
          Ncat(inNcat),
          omegamode(inomegamode),
          flataafitness(inflatfitness),
          flatcodonfitness(inflatcodonfitness),
          flatnucstat(inflatnucstat),
          flatnucrelrate(inflatnucrelrate) {
        init();
        Update();
    }

    ~AACodonMutSelMultipleOmega() override = default;

    void init() {
        blmode = 0;
        nucmode = 0;
        basemode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();
        Nstate = GetCodonStateSpace()->GetNstate();
        if (Ncat == -1) { Ncat = Nsite; }
        if (Ncat > Nsite) { Ncat = Nsite; }
        if (Ncat > 100) { Ncat = 100; }

        basemin = 0;
        if (baseNcat < 0) {
            basemin = 1;
            baseNcat = -baseNcat;
            if (baseNcat != 2) {
                std::cerr << "error in basencat\n";
                exit(1);
            }
        }

        if (flataafitness) {
            baseNcat = 1;
            Ncat = 1;
        }

        std::tuple<std::vector<std::vector<double>>, std::vector<size_t>> prefs{};
        if (!aaprofilesfile.empty() and aaprofilesfile != "Null") {
            if (flataafitness) {
                std::cerr << "Giving a preferences profiles file and setting 'flatfitness' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            prefs = aacodonmutsel::open_preferences(aaprofilesfile, Naa);
            clamp_profiles = true;
            Ncat = static_cast<int>(std::get<0>(prefs).size());
            if (static_cast<int>(std::get<1>(prefs).size()) == Nsite) {
                std::cout << "Input preferences profiles file has the same number of profiles "
                             "as the alignment, hence site allocations are considered fixed."
                          << std::endl;
                clamp_profiles_allocation = true;
            } else {
                std::cout << "Input preferences profiles file has a smaller number of profiles "
                             "than the alignment, hence site allocations are computed."
                          << std::endl;
            }
        }
        assert(Ncat <= Nsite);

        std::vector<double> opened_delta_omega_array;
        if (!delta_omega_array_file.empty() and delta_omega_array_file != "Null") {
            opened_delta_omega_array = aacodonmutsel::open_delta_omega_array(delta_omega_array_file);
            clamp_delta_omega_array = true;
            std::cout << "Input omega values are provided, hence are considered fixed but site "
                         "allocations are computed."
                      << std::endl;
            omegamode = 0;
            omegaNcat = opened_delta_omega_array.size();
        }
        std::cout << "Number of sites: " << Nsite << std::endl;
        std::cout << "Number of profile categories (cut-off for infinite mixture): " << Ncat
                  << '\n';
        std::cout << "Number of omega categories (finite mixture): " << omegaNcat << '\n';
        std::cout << "Additive shift in omega: " << omega_shift << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        Nbranch = tree->nb_branches();

        // Branch lengths
        blhypermean = 0.1;
        blhyperinvshape = 1.0;
        blhypermeanarray = new SimpleBranchArray<double>(*tree, blhypermean);
        branchlength = new GammaWhiteNoise(*tree, *blhypermeanarray, blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix: sampled from the prior, then overridden by the input
        // file (clamped) or set flat (clamped)
        std::vector<double> nucstat_{};
        if (!nucstatfile.empty() and nucstatfile != "Null") {
            if (flatnucstat) {
                std::cerr << "Giving a nucstat file and setting 'flatnucstat' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            nucstat_ = aacodonmutsel::open_nucstat(nucstatfile, Nnuc);
            clamp_nucstat = true;
        }
        std::vector<double> nucrelrate_{};
        if (!nucrelratefile.empty() and nucrelratefile != "Null") {
            if (flatnucrelrate) {
                std::cerr << "Giving a nucrelrate file and setting 'flatnucrelrate' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            nucrelrate_ = aacodonmutsel::open_nucrelrate(nucrelratefile, Nrr);
            clamp_nucrelrate = true;
        }

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, ((double)Nrr));
        if (clamp_nucrelrate) {
            std::copy(nucrelrate_.begin(), nucrelrate_.end(), nucrelrate.begin());
        } else if (flatnucrelrate) {
            nucrelrate.assign(Nrr, 1.0 / Nrr);
            clamp_nucrelrate = true;
        }

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, ((double)Nnuc));
        if (clamp_nucstat) {
            std::copy(nucstat_.begin(), nucstat_.end(), nucstat.begin());
        } else if (flatnucstat) {
            nucstat.assign(Nnuc, 1.0 / Nnuc);
            clamp_nucstat = true;
        }

        nucmatrix = new GTRSubMatrixCodon(Nnuc, nucrelrate, nucstat, true);

        // codonfitness
        codonfitness.assign(Nstate, 0);

        codonfitnesshypercenter.assign(Nstate, 1.0 / Nstate);
        codonfitnesshyperinvconc = 1.0 / Nstate;

        std::vector<double> codonprefs{};
        if (!codonfitnessfile.empty() and codonfitnessfile != "Null") {
            if (flatcodonfitness) {
                std::cerr
                    << "Giving a codon fitness file and setting 'flatcodonfitness' to true "
                       "is incompatible. These parameters are exclusive to one another."
                    << std::endl;
                exit(1);
            }
            codonprefs = aacodonmutsel::open_codonpreferences(codonfitnessfile, Nstate);
            clamp_codonfitness = true;
        }

        // synonymous groups from the genetic code (CodonStateSpace::Translation: codon -> AA)
        synonymousgroups.assign(Naa, std::vector<int>());
        for (int c = 0; c < Nstate; c++) {
            synonymousgroups[GetCodonStateSpace()->Translation(c)].push_back(c);
        }
        synonymousgroups.erase(
            std::remove_if(synonymousgroups.begin(), synonymousgroups.end(),
                [](const std::vector<int> &g) { return g.empty(); }),
            synonymousgroups.end());

        // prior draw: independent uniform Dirichlet(1,..,1) on each synonymous group
        // (single-codon AAs: theta = 1, no free parameter)
        for (auto const &g : synonymousgroups) {
            std::vector<double> sub(g.size(), 0.0);
            std::vector<double> center(g.size(), 1.0 / g.size());
            if (g.size() > 1) {
                Random::DirichletSample(sub, center, ((double)g.size()));
            } else {
                sub[0] = 1.0;
            }
            for (size_t k = 0; k < g.size(); k++) { codonfitness[g[k]] = sub[k]; }
        }
        if (clamp_codonfitness) {
            std::copy(codonprefs.begin(), codonprefs.end(), codonfitness.begin());
            NormaliseCodonFitness();  // file is a global simplex: rescale to sum 1 per AA
        } else if (flatcodonfitness) {
            // flat = uniform theta within each AA (theta_c = 1/n_a), i.e. codonweight = 1
            for (auto const &g : synonymousgroups) {
                for (int c : g) { codonfitness[c] = 1.0 / g.size(); }
            }
            clamp_codonfitness = true;
        }
        UpdateCodonWeights();

        // base distribution (can be skipped)
        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat, basekappa);
        baseprofile_occupancy = new OccupancySuffStat(baseNcat);

        basecenterhypercenter.assign(Naa, 1.0 / Naa);
        basecenterhyperinvconc = 1.0 / Naa;

        basecenterarray = new IIDDirichlet(baseNcat, basecenterhypercenter, basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Naa;
        baseconchyperinvshape = 1.0;
        baseconcentrationarray = new IIDGamma(baseNcat, baseconchypermean, baseconchyperinvshape);
        for (int k = 0; k < baseNcat; k++) { (*baseconcentrationarray)[k] = 20.0; }
        if (basemin == 1) { (*baseconcentrationarray)[0] = 1.0; }
        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat, Naa);
        componentalloc = new MultinomialAllocationVector(Ncat, baseweight->GetArray());
        componentcenterarray =
            new MixtureSelector<std::vector<double>>(basecenterarray, componentalloc);
        componentconcentrationarray =
            new MixtureSelector<double>(baseconcentrationarray, componentalloc);

        //
        // (truncated) Dirichlet mixture of aa fitness profiles
        //

        // Ncat fitness profiles iid from the base distribution
        componentaafitnessarray =
            new MultiDirichlet(componentcenterarray, componentconcentrationarray);
        if (flataafitness) { componentaafitnessarray->Flatten(); }

        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);

        // site allocations to the mixture (multinomial allocation)
        profile_alloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        // occupancy suff stats of site allocations for resampling weights (for both profile and
        // omega)
        profile_occupancy = new OccupancySuffStat(Ncat);
        omega_occupancy = new OccupancySuffStat(omegaNcat);

        // omega (fixed to 1 by default)
        delta_omegahypermean = 1.0;
        delta_omegahyperinvshape = 1.0;
        delta_omega_array = new IIDGamma(omegaNcat, delta_omegahypermean, delta_omegahyperinvshape);
        if (omegamode == 3 and omegaNcat == 1) { (*delta_omega_array)[0] = 1.0; }
        if (clamp_delta_omega_array) {
            for (int omegacat = 0; omegacat < omegaNcat; omegacat++) {
                (*delta_omega_array)[omegacat] = opened_delta_omega_array[omegacat];
            }
        }
        // will be a constant
        omega_weight_kappa = 1.0;
        omega_weight = new Dirichlet(omegaNcat, omega_weight_kappa);
        omega_alloc = new MultinomialAllocationVector(GetNsite(), omega_weight->GetArray());

        // selector, specifying which aa fitness array should be used for each site
        siteaafitnessarray =
            new MixtureSelector<std::vector<double>>(componentaafitnessarray, profile_alloc);

        if (clamp_profiles) {
            for (int cat = 0; cat < Ncat; cat++) {
                (*componentaafitnessarray)[cat] = std::get<0>(prefs)[cat];
            }
        }
        if (clamp_profiles_allocation) {
            for (int site = 0; site < Nsite; site++) {
                (*profile_alloc)[site] = static_cast<unsigned>(std::get<1>(prefs)[site]);
                assert(siteaafitnessarray->GetVal(site).size() == 20);
            }
        }

        // Ncat*omegaNcat mut sel codon matrices (based on the Ncat fitness profiles of the mixture,
        // and the omegaNcat finite mixture for omega)
        componentcodonmatrixbidimarray =
            new AACodonMutSelCodonMatrixBidimArray(GetCodonStateSpace(), nucmatrix, &codonweight,
                componentaafitnessarray, delta_omega_array, omega_shift);

        // Bidimselector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(
            componentcodonmatrixbidimarray, profile_alloc, omega_alloc);

        // Bidimselector, specifying which codon matrix should be used for each site (needed to
        // collect omegapathsuffstatarray)
        sitecodonsubmatrixarray = new DoubleMixtureSelector<AACodonMutSelOmegaCodonSubMatrix>(
            componentcodonmatrixbidimarray, profile_alloc, omega_alloc);

        phyloprocess =
            new PhyloProcess(tree.get(), codondata, branchlength, nullptr, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);

        siteomegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);
        componentomegapathsuffstatarray = new OmegaPathSuffStatArray(omegaNcat);
        tracer = std::unique_ptr<Tracer>(new Tracer(*this));
    }

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        if (blmode < 2) { model_node(info, "branchlength", *branchlength); }
        if (nucmode < 2) {
            model_node(info, "nucrelrate", nucrelrate);
            model_node(info, "nucstat", nucstat);
        }
        if (basemode < 2) {
            model_node(info, "basekappa", basekappa);
            model_node(info, "baseweight", *baseweight);
            model_node(info, "componentalloc", *componentalloc);
            model_node(info, "*basecenterarray", *basecenterarray);
            model_node(info, "*baseconcentrationarray", *baseconcentrationarray);
        }
        model_node(info, "kappa", kappa);
        model_node(info, "weight", *weight);
        model_node(info, "componentaafitnessarray", *componentaafitnessarray);
        model_node(info, "profile_alloc", *profile_alloc);
        model_node(info, "delta_omega_array", *delta_omega_array);
        model_node(info, "omega_alloc", *omega_alloc);
        model_node(info, "omega_weight", *omega_weight);
        model_node(info, "delta_omegahyperinvshape", delta_omegahyperinvshape);
        model_node(info, "delta_omegahypermean", delta_omegahypermean);
        model_node(info, "codonfitness", codonfitness);


        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        // 3x: per coding site (and not per nucleotide site)
        model_stat(info, "length", [this]() { return 3 * branchlength->GetTotalLength(); });
        model_stat(info, "ds", [this]() { return GetPredictedRelativedS(); });
        model_stat(info, "dnds", [this]() { return GetPredictedRelativedNdS(); });
        model_stat(info, "omega", [this]() { return GetMeanOmega(); });
        model_stat(info, "gc", [this]() { return GetGC(); });
        model_stat(info, "tstv", [this]() { return GetTsTv(); });
        model_stat(info, "nucsA", [this]() { return GetA(); });
        model_stat(info, "nucsC", [this]() { return GetC(); });
        model_stat(info, "nucsG", [this]() { return GetG(); });
        model_stat(info, "nucsT", [this]() { return GetT(); });
        model_stat(info, "nucrrAC", [this]() { return GetAC(); });
        model_stat(info, "nucrrAG", [this]() { return GetAG(); });
        model_stat(info, "nucrrAT", [this]() { return GetAT(); });
        model_stat(info, "nucrrCG", [this]() { return GetCG(); });
        model_stat(info, "nucrrCT", [this]() { return GetCT(); });
        model_stat(info, "nucrrGT", [this]() { return GetGT(); });
        model_stat(
            info, "omegaent", [this]() { return Random::GetEntropy(omega_weight->GetArray()); });
        model_stat(info, "ncluster", [this]() { return GetNcluster(); });
        model_stat(info, "kappa", kappa);
        if (baseNcat > 1) {
            if (basemin) {
                model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
                model_stat(info, "baseweight1", [this]() { return baseweight->GetVal(1); });
            } else {
                model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
                model_stat(info, "basekappa", basekappa);
            }
        }
        model_stat(info, "aaent", [this]() { return GetMeanAAEntropy(); });
        model_stat(info, "meanaaconc", [this]() { return GetMeanComponentAAConcentration(); });
        model_stat(info, "aacenterent", [this]() { return GetMeanComponentAAEntropy(); });
        model_stat(info, "statent", [this]() { return GetNucStatEntropy(); });
        model_stat(info, "rrent", [this]() { return GetNucRREntropy(); });
        model_stat(info, "codonent", [this]() { return GetCodonFitnessEntropy(); });
    }
    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    bool FlatFitness() const { return flataafitness; }

    bool FlatCodonFitness() const { return flatcodonfitness; }

    bool FlatNucStat() const { return flatnucstat; }

    bool FlatNucRelRate() const { return flatnucrelrate; }


    //! return branch length
    double GetBranchLength(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return branchlength->GetVal(tree->branch_index(node));
    }

    double GetNucRR(int i) {
        assert(i < Nrr);
        return nucrelrate[i];
    }
    double GetNucStat(int i) {
        assert(i < Nnuc);
        return nucstat[i];
    }

    //! relative preference theta_c among synonymous codons (sum 1 per AA)
    double GetCodonFitness(int i) {
        assert(i < Nstate);
        return codonfitness[i];
    }

    //! degeneracy-weighted preference n_a * theta_c used by the substitution
    //! matrices (SENCA); falls back to theta if not yet built.
    double GetCodonWeight(int i) {
        assert(i < Nstate);
        if (codonweight.size() == codonfitness.size()) {
            return codonweight[i];
        }
        return codonfitness[i];
    }

    double GetAASiteFitness(int site, int i) {
        assert(i < Naa);
        assert(site < Nsite);
        return siteaafitnessarray->GetVal(profile_alloc->GetVal(site))[i];
    }
    double GetAC() const { return nucrelrate[0]; }
    double GetAG() const { return nucrelrate[1]; }
    double GetAT() const { return nucrelrate[2]; }
    double GetCG() const { return nucrelrate[3]; }
    double GetCT() const { return nucrelrate[4]; }
    double GetGT() const { return nucrelrate[5]; }
    double GetA() const { return nucstat[0]; }
    double GetC() const { return nucstat[1]; }
    double GetG() const { return nucstat[2]; }
    double GetT() const { return nucstat[3]; }
    double GetGC() { return nucstat[1] + nucstat[2]; }
    double GetTsTv() {
        return (nucrelrate[1] + nucrelrate[4]) /
               (nucrelrate[0] + nucrelrate[2] + nucrelrate[3] + nucrelrate[5]);
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    int GetNtaxa() const { return taxonset->GetNtaxa(); }

    //! true if both models were built on the same alignment (same sites and taxa)
    bool IsDataEqual(AACodonMutSelMultipleOmega const &model) const {
        if (this->GetNsite() != model.GetNsite()) {
            std::cerr << "Non matching number of sites between the two models: " << this->GetNsite()
                      << " and " << model.GetNsite() << '\n';
            return false;
        }
        if (this->GetNtaxa() != model.GetNtaxa()) {
            std::cerr << "Non matching number of taxa between the two models: " << this->GetNtaxa()
                      << " and " << model.GetNtaxa() << '\n';
            return false;
        }
        for (int i{0}; i < GetNtaxa(); i++) {
            std::string taxa = taxonset->GetTaxon(i);
            if (not model.taxonset->TaxonPresent(taxa)) {
                std::cerr << "Taxon " << taxa << " is not present in the second model\n";
                return false;
            }
        }
        return true;
    }

    //! return current omega value for omega mixture of component k
    double GetComponentOmega(int k) const { return omega_shift + delta_omega_array->GetVal(k); }

    //! return current omega value for omega mixture of site
    double GetSiteOmega(int site) const { return GetComponentOmega(omega_alloc->GetVal(site)); }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void CorruptCodonMatrices() { componentcodonmatrixbidimarray->CorruptCodonMatrices(); }

    //! \brief tell codon matrices (of profile component k) that its parameters have changed and
    //! that it should be updated
    void CorruptProfileCodonMatrices(int k) {
        componentcodonmatrixbidimarray->CorruptRowCodonMatrices(k);
    }

    //! rescale codonfitness so that it sums to 1 within each synonymous group (SENCA codon layer)
    void NormaliseCodonFitness() {
        for (auto const &g : synonymousgroups) {
            double tot = 0;
            for (int c : g) { tot += codonfitness[c]; }
            for (int c : g) { codonfitness[c] /= tot; }
        }
    }

    //! codonweight_c = n_a * theta_c (degeneracy x intra-AA preference), used by the matrices
    void UpdateCodonWeights() {
        codonweight.resize(codonfitness.size());
        for (auto const &g : synonymousgroups) {
            for (int c : g) { codonweight[c] = g.size() * codonfitness[c]; }
        }
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateCodonWeights();
        UpdateNucMatrix();
        CorruptCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() {
        UpdateBaseOccupancies();
        UpdateProfileOccupancies();
        UpdateOmegaOccupancies();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(std::string name) {
        UpdateBaseOccupancies();
        UpdateOmegaOccupancies();
        UpdateProfileOccupancies();
        UpdateMatrices();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) { total += BranchLengthsLogPrior(); }
        if (nucmode < 2) { total += NucRatesLogPrior(); }
        total += CodonFitnessLogPrior();
        if (!clamp_profiles) {
            if (basemode < 2) {
                if (baseNcat > 1) {
                    total += BaseStickBreakingHyperLogPrior();
                    total += BaseStickBreakingLogPrior();
                }
                total += BaseLogPrior();
            }
            total += StickBreakingHyperLogPrior();
            total += StickBreakingLogPrior();
            total += AALogPrior();
        }
        if (omegamode < 2) { total += DeltaOmegaLogProb(); }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        return ret;
    }

    //! \brief log prior over delta_omega (for all components) + HyperLogPrior
    double DeltaOmegaLogProb() const { return DeltaOmegaLogPrior() + DeltaOmegaHyperLogPrior(); }

    //! \brief log prior over delta_omega (for all components)
    double DeltaOmegaLogPrior() const { return delta_omega_array->GetLogProb(); }

    //! \brief log prior over delta_omega (for component i)
    double DeltaOmegaLogPrior(int i) const { return delta_omega_array->GetLogProb(i); }

    //! \brief log prior over hyperparameter of prior over omega (here,
    //! delta_omegahyperinvshape ~ exponential of rate 10)
    //! delta_omegahypermean ~ exponential of rate 10)
    double DeltaOmegaHyperLogPrior() const {
        return -(delta_omegahyperinvshape + delta_omegahypermean) / 10;
    }

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over codonfitness: product over AA of uniform Dirichlet(1,..,1) on the n_a
    //! synonymous codons. Density w.r.t. Lebesgue measure on the (n_a - 1) free coordinates of
    //! each simplex = Gamma(n_a), hence log prior = sum_a logGamma(n_a) (constant; single-codon
    //! AAs contribute logGamma(1) = 0). Consistent with the Gamma representation used in
    //! MoveCodonGamma: theta = x / sum(x), x_c iid Gamma(1) per group.
    double CodonFitnessLogPrior() const {
        double total = 0;
        for (auto const &g : synonymousgroups) { total += Random::logGamma((double)g.size()); }
        return total;
    }

    //! log prior over concentration parameters basekappa of mixture of base
    //! distribution
    double BaseStickBreakingHyperLogPrior() const { return -basekappa / 10; }

    //! log prior over weights of stick breaking process of base distribution
    double BaseStickBreakingLogPrior() const { return baseweight->GetLogProb(basekappa); }

    //! log prior over concentration parameters kappa of stick-breaking mixture of
    //! amino-acid profiles
    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    //! log prior over weights of stick breaking process of amino-acid profiles
    double StickBreakingLogPrior() const { return weight->GetLogProb(kappa); }

    //! log prior over base center and concentration parameters
    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (std::isinf(total)) {
            std::cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    //! log prior over base center and concentration parameters of component k of
    //! base distribution
    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        return total;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const { return componentaafitnessarray->GetLogProb(); }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const { return componentaafitnessarray->GetLogProb(k); }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        double tot{0};
        for (int site{0}; site < GetNsite(); site++) {
            tot += sitepathsuffstatarray->GetVal(site).GetLogProb(sitesubmatrixarray->GetVal(site));
        }
        return tot;
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double PathSuffStatLogProb(int k) const {
        double tot{0};
        for (int site{0}; site < GetNsite(); site++) {
            if (profile_alloc->GetVal(site) == k) {
                tot += sitepathsuffstatarray->GetVal(site).GetLogProb(
                    sitesubmatrixarray->GetVal(site));
            }
        }
        return tot;
    }

    //! return log prob of the substitution mappings over sites allocated to omega
    //! component k of the omega mixture
    double PathSuffStatOmegaLogProb(int k) const {
        return componentomegapathsuffstatarray->GetVal(k).GetLogProb(GetComponentOmega(k));
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(
            basecenterarray->GetVal(k), 1.0 / baseconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving nucleotide mutation rate
    //! parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving codonfitness
    //! parameters (codon fitness)
    double CodonFitnessLogProb() const { return CodonFitnessLogPrior() + PathSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving aa hyper params (center and
    //! concentration) for component k of base distribution
    double BaseLogProb(int k) const { return BaseLogPrior(k) + BaseSuffStatLogProb(k); }

    //! log prob factor to be recomputed when moving basekappa, the concentration
    //! parameter of the second-level strick breaking process (base distribution)
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
    }

    //! log prob factor to be recomputed when moving kappa, the concentration
    //! parameter of the first-level strick breaking process
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics of substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! collect omega sufficient statistics of substitution mappings across sites
    void CollectSiteOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonsubmatrixarray, *sitepathsuffstatarray);
    }

    //! gather site-specific omega sufficient statistics component-wise (per omega component)
    void CollectComponentOmegaPathSuffStat() {
        componentomegapathsuffstatarray->Clear();
        componentomegapathsuffstatarray->Add(*siteomegapathsuffstatarray, *omega_alloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //-------------------
    //  Moves
    //-------------------

    //! complete MCMC move schedule
    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        accept_iter++;
        AcceptanceTracker::get().periodicReport(accept_iter);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            if (blmode < 2) { MoveBranchLengths(); }

            CollectSitePathSuffStat();

            // clamped parameters (loaded from a file, or flat) are never moved
            if (nucmode < 2) {
                if (!clamp_nucrelrate) { MoveNucRR(); }
                if (!clamp_nucstat and !clamp_codonfitness) {
                    MoveCodonNucStat();
                    MoveCodonFull();
                } else if (!clamp_nucstat) {
                    MoveNucStat();
                } else if (!clamp_codonfitness) {
                    MoveCodonFull();
                }
            }

            if (omegamode < 2) { MoveOmegaMixture(3); }

            if (!flataafitness) {
                aachrono.Start();
                MoveAAMixture(3);
                aachrono.Stop();

                basechrono.Start();
                if (basemode < 2 and !clamp_profiles) { MoveBase(3); }
                basechrono.Stop();
            }
            totchrono.Stop();
        }
    }

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1) { ResampleBaseAlloc(); }
        MoveBaseMixture(nrep);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() { ResampleBranchLengths(); }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        MoveNucRR();
        MoveNucStat();
    }

    void MoveNucRR() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &AACodonMutSelMultipleOmega::NucRatesLogProb,
            &AACodonMutSelMultipleOmega::UpdateMatrices, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &AACodonMutSelMultipleOmega::NucRatesLogProb,
            &AACodonMutSelMultipleOmega::UpdateMatrices, this);
    }

    void MoveNucStat() {
        Move::Profile(nucstat, 0.1, 1, 2, &AACodonMutSelMultipleOmega::NucRatesLogProb,
            &AACodonMutSelMultipleOmega::UpdateMatrices, this);
        Move::Profile(nucstat, 0.01, 3, 2, &AACodonMutSelMultipleOmega::NucRatesLogProb,
            &AACodonMutSelMultipleOmega::UpdateMatrices, this);
    }

    //! MH move on codonfitnesses and nucstat parameters
    double MoveCodonNucStat(
        double tuning_codon, int n_codon, double tuning_stat, int n_stat, int nrep) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> BKCodonFitness(codonfitness.size(), 0);
        std::vector<double> BKNucStat(nucstat.size(), 0);
        for (int rep = 0; rep < nrep; rep++) {
            BKCodonFitness = codonfitness;
            BKNucStat = nucstat;

            // the path log prob is common to both parameters: count it once
            double deltalogprob =
                -CodonFitnessLogPrior() - NucRatesLogPrior() - PathSuffStatLogProb();
            deltalogprob += ProposeSynonymousPairMove(tuning_codon, n_codon);
            deltalogprob += Random::ProfileProposeMove(nucstat, nucstat.size(), tuning_stat, n_stat);
            UpdateMatrices();
            deltalogprob += CodonFitnessLogPrior() + NucRatesLogPrior() + PathSuffStatLogProb();

            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                codonfitness = BKCodonFitness;
                nucstat = BKNucStat;
                UpdateMatrices();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    void MoveCodonNucStat() {
        int nrep{3};
        MoveCodonNucStat(1, 1, 0.1, 1, nrep);
        MoveCodonNucStat(0.1, 3, 0.1, 1, nrep);
        MoveCodonNucStat(0.01, 3, 0.01, 1, nrep);
        MoveCodonNucStat(0.001, 3, 0.01, 1, nrep);
    }
    //! "vase communiquants" proposal restricted to synonymous codons: n times, draw a
    //! synonymous group with >= 2 codons uniformly, then an ordered pair (i1,i2) of distinct
    //! codons in it, and redistribute theta[i1]+theta[i2] between them (reflected uniform
    //! slide, as Random::ProfileProposeMove). The pair sum is conserved, so sum_{syn(a)} theta
    //! = 1 holds exactly for every AA (the global ProfileProposeMove over the 61 entries would
    //! only conserve the global sum and move mass between AAs). The kernel is symmetric
    //! (group choice and ordered-pair choice do not depend on the state): log Hastings = 0.
    double ProposeSynonymousPairMove(double tuning, int n) {
        std::vector<int> multi;  // groups with at least 2 codons
        for (size_t g = 0; g < synonymousgroups.size(); g++) {
            if (synonymousgroups[g].size() > 1) { multi.push_back(g); }
        }
        for (int k = 0; k < n; k++) {
            auto const &grp = synonymousgroups[multi[(int)(Random::Uniform() * multi.size())]];
            int m = grp.size();
            int a = (int)(Random::Uniform() * m);
            int b = (int)(Random::Uniform() * (m - 1));
            if (b >= a) { b++; }
            int i1 = grp[a];
            int i2 = grp[b];
            double tot = codonfitness[i1] + codonfitness[i2];
            double x = codonfitness[i1] + tot * tuning * (Random::Uniform() - 0.5);
            while ((x < 0) || (x > tot)) {
                if (x < 0) { x = -x; }
                if (x > tot) { x = 2 * tot - x; }
            }
            codonfitness[i1] = x > 1e-50 ? x : 1e-50;
            codonfitness[i2] = tot - x > 1e-50 ? tot - x : 1e-50;
        }
        return 0;
    }

    //! MH move on codonfitnesses parameters (group-wise additive pair moves)
    void MoveCodon() {
        MoveCodonAddPair(1.0, 1, 3);
        MoveCodonAddPair(0.1, 3, 3);
        MoveCodonAddPair(0.01, 3, 3);
        MoveCodonAddPair(0.001, 3, 3);
    }

    //! MH move on codonfitness: additive compensated move on n synonymous pairs of entries
    double MoveCodonAddPair(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(codonfitness.size(), 0);
        for (int rep = 0; rep < nrep; rep++) {
            bk = codonfitness;
            double deltalogprob = -CodonFitnessLogPrior() - PathSuffStatLogProb();
            deltalogprob += ProposeSynonymousPairMove(tuning, n);
            UpdateMatrices();
            deltalogprob += CodonFitnessLogPrior() + PathSuffStatLogProb();

            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                codonfitness = bk;
                UpdateMatrices();
            }
            ntot++;
        }
        AcceptanceTracker::get().add("aac_pair_t" + std::to_string((int)(tuning * 1000.0)) + "_n" + std::to_string(n), nacc, ntot);
        return nacc / ntot;
    }

    //! schedule of additive pair moves (as CompMoveCodonProfiles in mutselc)
    void MoveCodonAddPairSchedule() {
        MoveCodonAddPair(0.3, 2, 1);
        MoveCodonAddPair(0.15, 4, 1);
        MoveCodonAddPair(0.07, 8, 1);
        MoveCodonAddPair(0.03, 16, 1);
        MoveCodonAddPair(0.01, 32, 1);
        MoveCodonAddPair(1, 2, 1);
    }

    //! helper function: log density of the Gamma(1) variables x_c underlying each
    //! synonymous Dirichlet(1,..,1) (alpha = 1: log density = -x_c). Single-codon AAs have no
    //! free parameter and are excluded.
    double GammaCodonLogPrior(const std::vector<double> &x) {
        double total = 0;
        for (auto const &g : synonymousgroups) {
            if (g.size() > 1) {
                for (int c : g) { total -= x[c]; }
            }
        }
        return total;
    }

    //! MH move on codonfitness: multiplicative move, per synonymous group (Gamma
    //! representation of each Dirichlet: theta_c = x_c / sum_{syn} x, x_c ~ Gamma(1)).
    //! Dividing by the group sum after the move keeps sum_{syn(a)} theta = 1 for every AA.
    //! As in mutselc, the auxiliary sum z_g ~ Gamma(n_a) is redrawn (Gibbs) so x = z_g * theta.
    double MoveCodonGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;

        std::vector<double> x(Nstate, 0);
        for (auto const &g : synonymousgroups) {
            double z = Random::sGamma((double)g.size());
            for (int c : g) { x[c] = z * codonfitness[c]; }
        }

        std::vector<double> bkx = x;
        std::vector<double> bkcodon = codonfitness;

        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -GammaCodonLogPrior(x) - PathSuffStatLogProb();

            double loghastings = 0;
            for (auto const &g : synonymousgroups) {
                if (g.size() == 1) { continue; }  // theta = 1 fixed
                double z = 0;
                for (int c : g) {
                    double m = tuning * (Random::Uniform() - 0.5);
                    x[c] *= exp(m);
                    z += x[c];
                    loghastings += m;
                }
                for (int c : g) {
                    codonfitness[c] = x[c] / z;
                    if (codonfitness[c] < 1e-50) { codonfitness[c] = 1e-50; }
                }
            }

            deltalogprob += loghastings;

            UpdateMatrices();

            deltalogprob += GammaCodonLogPrior(x) + PathSuffStatLogProb();

            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
                bkcodon = codonfitness;
                bkx = x;
            } else {
                codonfitness = bkcodon;
                x = bkx;
                UpdateMatrices();
            }
            ntot++;
        }
        AcceptanceTracker::get().add("aac_codongamma", nacc, ntot);
        return nacc / ntot;
    }

    //! full move schedule on codonfitness: additive pairs then multiplicative
    void MoveCodonFull() {
        MoveCodonAddPairSchedule();
        MoveCodonGamma(0.1, 1);
        MoveCodonGamma(0.1, 1);
    }

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            if (!clamp_profiles) {
                MoveAAProfiles();
                ResampleEmptyProfileComponents();
            }
            if (!clamp_profiles_allocation) {
                ResampleProfileAlloc();
                ProfileLabelSwitchingMove();
                ResampleProfileWeights();
                MoveKappa();
            }
            if (!clamp_profiles_allocation or !clamp_profiles) { CorruptCodonMatrices(); }
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyProfileComponents() {
        componentaafitnessarray->PriorResample(*profile_occupancy);
        componentcodonmatrixbidimarray->CorruptCodonMatricesRowOccupancy(*profile_occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        MoveAA(1.0, 1, nrep);
        MoveAA(0.1, 3, nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        MoveAAGamma(3.0, nrep);
        MoveAAGamma(1.0, nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveAA(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                std::vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) { bk[l] = aa[l]; }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    CorruptProfileCodonMatrices(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Naa; l++) { aa[l] = bk[l]; }
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(
        const std::vector<double> &x, const std::vector<double> &aacenter, double aaconc) {
        double total = 0;
        for (int l = 0; l < Naa; l++) {
            total += (aaconc * aacenter[l] - 1) * log(x[l]) - x[l] -
                     Random::logGamma(aaconc * aacenter[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MoveAAGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                double aaconc = componentconcentrationarray->GetVal(i);
                const std::vector<double> &aacenter = componentcenterarray->GetVal(i);

                std::vector<double> &aa = (*componentaafitnessarray)[i];
                std::vector<double> x(Naa, 0);
                double z = Random::sGamma(aaconc);
                for (int l = 0; l < Naa; l++) { x[l] = z * aa[l]; }

                double bkz = z;
                std::vector<double> bkx = x;
                std::vector<double> bkaa = aa;

                for (int rep = 0; rep < nrep; rep++) {
                    double deltalogprob =
                        -GammaAALogPrior(x, aacenter, aaconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l = 0; l < Naa; l++) {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l = 0; l < Naa; l++) {
                        aa[l] = x[l] / z;
                        if (aa[l] < 1e-50) { aa[l] = 1e-50; }
                    }

                    deltalogprob += loghastings;

                    CorruptProfileCodonMatrices(i);

                    deltalogprob += GammaAALogPrior(x, aacenter, aaconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    } else {
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! Gibbs resample mixture allocations
    void ResampleProfileAlloc() {
        std::vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetProfileAllocPostProb(i, postprob);
            profile_alloc->GibbsResample(i, postprob);
        }
        UpdateProfileOccupancies();
    }

    //! update mixture profile occupancy suff stats (for resampling mixture weights)
    void UpdateProfileOccupancies() {
        profile_occupancy->Clear();
        profile_occupancy->AddSuffStat(*profile_alloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetProfileAllocPostProb(int site, std::vector<double> &postprob) {
        double max = 0;
        const std::vector<double> &w = weight->GetArray();
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        for (int i = 0; i < Ncat; i++) {
            double tmp = suffstat.GetLogProb(
                componentcodonmatrixbidimarray->GetVal(i, omega_alloc->GetVal(site)));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Ncat; i++) { postprob[i] /= total; }
    }

    //! MCMC sequence for label switching moves
    void ProfileLabelSwitchingMove() {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5, *profile_occupancy, permut);
        profile_alloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    //! Gibbs resample mixture weights (based on profile_occupancy suff stats)
    void ResampleProfileWeights() { weight->GibbsResample(*profile_occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        Move::Scaling(kappa, 1.0, 10, &AACodonMutSelMultipleOmega::StickBreakingHyperLogProb,
            &AACodonMutSelMultipleOmega::NoUpdate, this);
        Move::Scaling(kappa, 0.3, 10, &AACodonMutSelMultipleOmega::StickBreakingHyperLogProb,
            &AACodonMutSelMultipleOmega::NoUpdate, this);
    }

    //! MCMC module for the base mixture
    void MoveBaseMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseComponents(10);
            ResampleBaseEmptyComponents();
            if (baseNcat > 1) {
                if (!basemin) { BaseLabelSwitchingMove(); }
                ResampleBaseWeights();
                MoveBaseKappa();
            }
        }
    }

    //! MCMC module for moving the center and concentration parameters of the
    //! components of the the base mixture
    void MoveBaseComponents(int nrep) {
        CollectBaseSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseCenters(1.0, 1);
            MoveBaseCenters(1.0, 3);
            MoveBaseCenters(0.3, 3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.3);
        }
    }

    //! collect suff stats for moving center and concentration parameters of the
    //! base mixture
    void CollectBaseSuffStat() {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray, *componentalloc);
    }

    //! MCMC module for moving the center parameters of the components of the the
    //! base mixture
    double MoveBaseCenters(double tuning, int n) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(Naa, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseprofile_occupancy->GetVal(k)) {
                std::vector<double> &aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    aa = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MCMC module for moving the concentration parameters of the components of
    //! the the base mixture
    double MoveBaseConcentrations(double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int k = basemin; k < baseNcat; k++) {
            if (baseprofile_occupancy->GetVal(k)) {
                double &c = (*baseconcentrationarray)[k];
                double bk = c;
                double deltalogprob = -BaseLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                c *= e;
                deltalogprob += m;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    c = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! resample empty components of the base mixture from the prior
    void ResampleBaseEmptyComponents() {
        basecenterarray->PriorResample(*baseprofile_occupancy);
        baseconcentrationarray->PriorResample(*baseprofile_occupancy);
        if (basemin == 1) { (*baseconcentrationarray)[0] = 1.0; }
    }

    //! Gibbs resample base mixture allocations
    void ResampleBaseAlloc() {
        std::vector<double> postprob(baseNcat, 0);
        for (int i = 0; i < Ncat; i++) {
            GetBaseAllocPostProb(i, postprob);
            componentalloc->GibbsResample(i, postprob);
            if ((componentalloc->GetVal(i) < 0) || (componentalloc->GetVal(i) >= baseNcat)) {
                std::cerr << "error in ResampleBaseAlloc: out of bound\n";
                exit(1);
            }
        }
        UpdateBaseOccupancies();
    }

    //! update base profile_occupancy suff stats (for moving base weights)
    void UpdateBaseOccupancies() {
        baseprofile_occupancy->Clear();
        baseprofile_occupancy->AddSuffStat(*componentalloc);
    }

    //! get allocation posterior probability of a component of the first-level
    //! mixture to the components of the second-level mixture
    void GetBaseAllocPostProb(int cat, std::vector<double> &postprob) {
        double max = 0;
        const std::vector<double> &w = baseweight->GetArray();
        for (int i = 0; i < baseNcat; i++) {
            double tmp = Random::logDirichletDensity(componentaafitnessarray->GetVal(cat),
                basecenterarray->GetVal(i), baseconcentrationarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < baseNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < baseNcat; i++) { postprob[i] /= total; }
    }

    //! MCMC sequence for label switching moves of the base mixture
    void BaseLabelSwitchingMove() {
        Permutation permut(baseNcat);
        baseweight->LabelSwitchingMove(5, *baseprofile_occupancy, permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    //! Gibbs resample base mixture weights (based on profile_occupancy suff stats)
    void ResampleBaseWeights() { baseweight->GibbsResample(*baseprofile_occupancy); }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa() {
        Move::Scaling(basekappa, 1.0, 10,
            &AACodonMutSelMultipleOmega::BaseStickBreakingHyperLogProb,
            &AACodonMutSelMultipleOmega::NoUpdate, this);
        Move::Scaling(basekappa, 0.3, 10,
            &AACodonMutSelMultipleOmega::BaseStickBreakingHyperLogProb,
            &AACodonMutSelMultipleOmega::NoUpdate, this);
    }

    //! MCMC module for the mixture of omega
    void MoveOmegaMixture(int nrep) {
        CollectSiteOmegaPathSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            ResampleOmegaAlloc();
            UpdateOmegaOccupancies();
            if (!clamp_delta_omega_array) {
                CollectComponentOmegaPathSuffStat();
                MoveOmegaValues(1.0, 3);
                MoveOmegaValues(0.1, 3);
                ResampleEmptyOmegaComponents();
                MoveOmegaHyper();
            }
            ResampleOmegaWeights();
        }
        CorruptCodonMatrices();
    }

    //! MH move on omega (per component): scaling move
    double MoveOmegaValues(double tuning, int nrep) {
        double bk{0}, nacc{0}, ntot{0};
        for (int i = 0; i < omegaNcat; i++) {
            if (omega_occupancy->GetVal(i)) {
                for (int rep = 0; rep < nrep; rep++) {
                    bk = delta_omega_array->GetVal(i);

                    double deltalogprob = -(PathSuffStatOmegaLogProb(i) + DeltaOmegaLogPrior(i));
                    double loghastings = tuning * (Random::Uniform() - 0.5);
                    double hastings = exp(loghastings);

                    (*delta_omega_array)[i] *= hastings;

                    deltalogprob += loghastings;
                    deltalogprob += PathSuffStatOmegaLogProb(i) + DeltaOmegaLogPrior(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        (*delta_omega_array)[i] = bk;
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    void ResampleEmptyOmegaComponents() { delta_omega_array->PriorResample(*omega_occupancy); };

    //! update mixture omega occupancy suff stats (for resampling mixture weights)
    void UpdateOmegaOccupancies() {
        omega_occupancy->Clear();
        omega_occupancy->AddSuffStat(*omega_alloc);
    }

    //! Gibbs resample omega mixture allocations
    void ResampleOmegaAlloc() {
        std::vector<double> postprob(omegaNcat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetOmegaAllocPostProb(i, postprob);
            omega_alloc->GibbsResample(i, postprob);
        }
    }

    //! get omega allocation posterior probabilities for a given site
    double GetOmegaAllocPostProb(int site, std::vector<double> &postprob) {
        double max = 0;

        const std::vector<double> &w = omega_weight->GetArray();
        const OmegaPathSuffStat &suffstat = siteomegapathsuffstatarray->GetVal(site);

        for (int i = 0; i < omegaNcat; i++) {
            double tmp = suffstat.GetLogProb(GetComponentOmega(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) { max = tmp; }
        }

        double total = 0;
        for (int i = 0; i < omegaNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < omegaNcat; i++) { postprob[i] /= total; }

        return log(total) + max;
    }

    //! Gibbs resample omega mixture weights (based on omega occupancy suff stats)
    void ResampleOmegaWeights() { omega_weight->GibbsResample(*omega_occupancy); }


    //! Scaling moves on omega hyperparameters
    void MoveOmegaHyper() {
        if (omegaNcat > 1) {
            Move::Scaling(delta_omegahypermean, 1.0, 10,
                &AACodonMutSelMultipleOmega::DeltaOmegaLogProb,
                &AACodonMutSelMultipleOmega::NoUpdate, this);
            Move::Scaling(delta_omegahypermean, 0.3, 10,
                &AACodonMutSelMultipleOmega::DeltaOmegaLogProb,
                &AACodonMutSelMultipleOmega::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 1.0, 10,
                &AACodonMutSelMultipleOmega::DeltaOmegaLogProb,
                &AACodonMutSelMultipleOmega::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 0.3, 10,
                &AACodonMutSelMultipleOmega::DeltaOmegaLogProb,
                &AACodonMutSelMultipleOmega::NoUpdate, this);
        }
    }


    //-------------------
    // Traces and Monitors
    // ------------------

    //! return mutation rate between nucleotides from the mutation matrix
    double GetNucRate(int i, int j) const {
        nucmatrix->UpdateMatrix();
        // full Q entry (nucrate * pi_j), not the symmetric exchangeability: matches the
        // --nuc output of the other mut-sel readers (see mutselc fix d57d4b4)
        return (*nucmatrix)(i, j);
    }


    //! return number of occupied components in first-level mixture (mixture of
    //! amino-acid fitness profiles)
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) { n++; }
        }
        return n;
    }

    //! return number of occupied components in base distribution
    int GetBaseNcluster() const {
        int n = 0;
        for (int i = 0; i < baseNcat; i++) {
            if (baseprofile_occupancy->GetVal(i)) { n++; }
        }
        return n;
    }


    //! return mean omega
    double GetMeanOmega() const {
        double tot = 0;
        for (int i{0}; i < omegaNcat; i++) { tot += GetComponentOmega(i); }
        return tot / omegaNcat;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanAAEntropy() const { return componentaafitnessarray->GetMeanEntropy(); }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentAAConcentration() const {
        double tot = 0;
        double totw = 0;
        for (int i = basemin; i < baseNcat; i++) {
            tot += baseprofile_occupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
            totw += baseprofile_occupancy->GetVal(i);
        }
        return tot / totw;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentAAEntropy() const {
        double tot = 0;
        for (int i = 0; i < baseNcat; i++) {
            tot +=
                baseprofile_occupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucstat); }

    //! return entropy of vector of codon fitnesses
    double GetCodonFitnessEntropy() const { return Random::GetEntropy(codonfitness); }


    //! return predicted omega induced by the mutation-selection balance
    double GetPredictedOmegaKnot() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) { mean += GetPredictedSiteOmegaKnot(i); }
        mean /= GetNsite();
        return mean;
    }

    //! return predicted omega at a given site induced by the mutation-selection balance
    double GetPredictedSiteOmegaKnot(int site) const {
        return sitecodonsubmatrixarray->GetVal(site).GetPredictedDNDS();
    }

    //! return effective dnds taking into the mutation-selection balance and omega
    double GetPredictedEffectivedNdS() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += GetSiteOmega(i) * sitecodonsubmatrixarray->GetVal(i).GetPredictedDNDS();
        }
        mean /= GetNsite();
        return mean;
    }

    double GetPredictedEffectivedS() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += sitecodonsubmatrixarray->GetVal(i).GetPredictedDS();
        }
        mean /= GetNsite();
        return mean;
    }

    double GetPredictedRelativedNdS() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += GetSiteOmega(i) * sitecodonsubmatrixarray->GetVal(i).GetPredictedRelativeDNDS();
        }
        mean /= GetNsite();
        return mean;
    }

    double GetPredictedRelativedN() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += GetSiteOmega(i) * sitecodonsubmatrixarray->GetVal(i).GetPredictedRelativeDN();
        }
        mean /= GetNsite();
        return mean;
    }

    double GetPredictedRelativedS() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += sitecodonsubmatrixarray->GetVal(i).GetPredictedRelativeDS();
        }
        mean /= GetNsite();
        return mean;
    }

    const Tree &GetTree() const { return *tree; }

    const int GetNcat() const { return Ncat; }

    const int GetOmegaNcat() const { return omegaNcat; }

    const double GetProfileAA(int i, int a) const { return componentaafitnessarray->GetVal(i)[a]; }
    const int GetProfileAlloc(int i) const { return profile_alloc->GetVal(i); }
    const double GetOmega(int i) const { return delta_omega_array->GetVal(i); }
    const int GetOmegaAlloc(int i) const { return omega_alloc->GetVal(i); }

    const std::vector<double> &GetProfile(int i) const { return siteaafitnessarray->GetVal(i); }

    static std::string GetModelName() { return "AACodonMutSelMultipleOmega"; }

    //! stamp without the nucstat/nucrelrate file names (used for the .pvalues output)
    void GetPvaluesStamp(std::ostream &os) const {
        os << GetModelName() << '\t';
        os << datafile << '\t' << treefile << '\t' << aaprofilesfile << '\t' << codonfitnessfile
           << '\t';
        os << delta_omega_array_file << '\t';
        os << omegaNcat << '\t';
        {
            // omega_shift is the only floating-point value of the stamp: write it with full
            // precision so that it survives a resume exactly; restore the stream state since
            // the tracer line is written on the same stream afterwards.
            auto old_precision = os.precision(std::numeric_limits<double>::max_digits10);
            os << omega_shift << '\t';
            os.precision(old_precision);
        }
        // init() makes baseNcat positive: the sign carries basemin through a resume
        os << (basemin ? -baseNcat : baseNcat) << '\t';
        os << Ncat << '\t';
        os << omegamode << '\t';
        os << flataafitness << '\t';
        os << flatcodonfitness << '\t';
        os << flatnucstat << '\t';
        os << flatnucrelrate << '\t';
    }

    //! full stamp, written to the .param file for resume: the pvalues stamp, followed by the
    //! nucstat/nucrelrate file names (without them, a resumed chain would not know that these
    //! parameters are clamped)
    void GetModelStamp(std::ostream &os) const {
        GetPvaluesStamp(os);
        os << nucstatfile << '\t' << nucrelratefile << '\t';
    }

    void ToStream(std::ostream &os) const {
        GetModelStamp(os);
        tracer->write_line(os);
    }

    AACodonMutSelMultipleOmega(std::istream &is) {
        // the stamp is tab-delimited (see GetPvaluesStamp), so that file paths may contain spaces
        auto next_field = [&is]() {
            std::string field;
            if (!std::getline(is >> std::ws, field, '\t')) {
                std::cerr << "error: truncated model stamp in checkpoint\n";
                exit(1);
            }
            return field;
        };
        auto next_number = [&next_field](double &value) { value = std::stod(next_field()); };
        std::string model_name = next_field();
        if (model_name != GetModelName()) {
            std::cerr << "Expected " << GetModelName() << " for model name, got " << model_name
                      << "\n";
            exit(1);
        }
        datafile = next_field();
        treefile = next_field();
        aaprofilesfile = next_field();
        codonfitnessfile = next_field();
        delta_omega_array_file = next_field();
        double value{0};
        next_number(value);
        omegaNcat = static_cast<int>(value);
        next_number(omega_shift);
        next_number(value);
        baseNcat = static_cast<int>(value);
        next_number(value);
        Ncat = static_cast<int>(value);
        next_number(value);
        omegamode = static_cast<int>(value);
        next_number(value);
        flataafitness = value != 0;
        next_number(value);
        flatcodonfitness = value != 0;
        next_number(value);
        flatnucstat = value != 0;
        next_number(value);
        flatnucrelrate = value != 0;
        nucstatfile = next_field();
        nucrelratefile = next_field();
        init();
        tracer->read_line(is);
        Update();
    };
};
