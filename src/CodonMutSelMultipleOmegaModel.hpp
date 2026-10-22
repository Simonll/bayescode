#pragma once

#include "Chrono.hpp"
#include "CodonMutSelCodonMatrixBidimArray.hpp"
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

#define CLASS_NAME(class_name) #class_name

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- double Dirichlet process version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific codon fitness profiles F_ic, for site i and
 * codon c
 * - an array of site-specific omega multiplier, capturing deviations of the non-syn rate from the
 * model; this parameter is fixed to 1 by default.
 *
 * Site-specific codon fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration
 * parameter kappa, and truncated at Ncat.
 *
 * The base distribution of this Dirichlet process, G_0, is itself a (truncated
 * stick breaking) mixture of Dirichlet distributions, parameterized by
 * basekappa, and truncated at baseNcat. Each component of this mixture
 * Dirichlet is parameterized by a center (a 61-freqeuncy vector) and a
 * concentration.
 *
 * Concerning the first-level stick-breaking process, by default, Ncat == 100
 * (can be changed using the -ncat option). As for baseNcat, it is equal to 1,
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
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 61.
 *
 */


//! Open a file
inline std::vector<double> open_delta_omega_array(std::string const &file_name) {
    std::vector<double> delta_omega_array{};
    std::ifstream input_stream(file_name);
    if (!input_stream) {
        std::cerr << "Delta-omega array file " << file_name << " doesn't exist" << std::endl;
        exit(1);
    }

    std::string line;
    while (getline(input_stream, line)) {
        if (!line.empty()) { delta_omega_array.push_back(stod(line)); }
    }
    return delta_omega_array;
}

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

inline std::vector<double> open_nucstat(std::string const &file_name, int size) {
    return open_simplex_vector(file_name, size, "nucstat");
}

inline std::vector<double> open_nucrelrate(std::string const &file_name, int size) {
    return open_simplex_vector(file_name, size, "nucrelrate");
}
class CodonMutSelMultipleOmegaModel : public ChainComponent {
    std::string datafile, treefile, profilesfile{"Null"}, nucstatfile{"Null"},
        nucrelratefile{"Null"};
    bool clamp_profiles{false}, clamp_profiles_allocation{false}, clamp_nucstat{false},
        clamp_nucrelrate{false}, loaded_nucstat{false}, loaded_nucrelrate{false},
        loaded_profiles{false}, loaded_delta_omega_array{false}, loaded_profiles_allocation{false};
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

    MultiDirichlet *componentcodonfitnessarray;
    DirichletSuffStatArray *basesuffstatarray;

    MultinomialAllocationVector *profile_alloc;
    // sites allocated to each profile component; kept in sync with profile_alloc by
    // UpdateProfileSites() so that per-component log probs cost O(sites of k), not O(Nsite)
    std::vector<std::vector<int>> profile_sites;

    // an array of codon matrices (one for each distinct aa fitness profile)
    CodonMutSelCodonMatrixBidimArray *componentcodonmatrixbidimarray;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;
    DoubleMixtureSelector<CodonMutSelOmegaCodonSubMatrix> *sitecodonsubmatrixarray;

    MixtureSelector<std::vector<double>> *sitecodonfitnessarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;

    // 0: free without shrinkage
    // 2: shared across genes
    // 3: fixed
    // (modes 1 and 2 are not implemented for this model)

    // currently: free without shrinkage
    int blmode;
    // currently: free without shrinkage
    int nucmode;
    // currently: free without shrinkage
    int basemode;

    // currently: free without shrinkage (0) or fixed (3)
    int omegamode;

    bool flatfitness;
    bool flatnucstat;
    bool flatnucrelrate;
    Chrono profilechrono;
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
    CodonMutSelMultipleOmegaModel(std::string indatafile, std::string intreefile,
        std::string inprofilesfile, std::string innucstat, std::string innucrelrate,
        int inomegamode, int inNcat, int inbaseNcat, int inomegaNcat, double inomegashift,
        bool inflatfitness, std::string indelta_omega_array_file, bool inflatnucstat,
        bool inflatnucrelrate)
        : datafile(std::move(indatafile)),
          treefile(std::move(intreefile)),
          profilesfile(std::move(inprofilesfile)),
          nucstatfile(std::move(innucstat)),
          nucrelratefile(std::move(innucrelrate)),
          delta_omega_array_file(std::move(indelta_omega_array_file)),
          omegaNcat(inomegaNcat),
          omega_shift(inomegashift),
          baseNcat(inbaseNcat),
          Ncat(inNcat),
          omegamode(inomegamode),
          flatfitness(inflatfitness),
          flatnucstat(inflatnucstat),
          flatnucrelrate(inflatnucrelrate) {
        init();
        Update();
    }

    ~CodonMutSelMultipleOmegaModel() override {
        // delete in reverse order of creation in init(); member destructors do not
        // delete their cross-references, so this order is safe.
        delete componentomegapathsuffstatarray;
        delete siteomegapathsuffstatarray;
        delete sitepathsuffstatarray;
        delete phyloprocess;
        delete sitecodonsubmatrixarray;
        delete sitesubmatrixarray;
        delete componentcodonmatrixbidimarray;
        delete sitecodonfitnessarray;
        delete omega_alloc;
        delete omega_weight;
        delete delta_omega_array;
        delete omega_occupancy;
        delete profile_occupancy;
        delete profile_alloc;
        delete weight;
        delete componentcodonfitnessarray;
        delete componentconcentrationarray;
        delete componentcenterarray;
        delete componentalloc;
        delete basesuffstatarray;
        delete baseconcentrationarray;
        delete basecenterarray;
        delete baseprofile_occupancy;
        delete baseweight;
        delete nucmatrix;
        delete lengthpathsuffstatarray;
        delete branchlength;
        delete blhypermeanarray;
        delete codondata;
        delete data;
    }

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
        if (Ncat > 100) {
            std::cerr << "warning: ncat capped at 100 (requested " << Ncat << ")\n";
            Ncat = 100;
        }

        basemin = 0;
        if (baseNcat < 0) {
            basemin = 1;
            baseNcat = -baseNcat;
            if (baseNcat != 2) {
                std::cerr << "error in basencat\n";
                exit(1);
            }
        }

        if (flatfitness) {
            baseNcat = 1;
            basemin = 0;
            Ncat = 1;
        }

        std::tuple<std::vector<std::vector<double>>, std::vector<size_t>> prefs{};
        if (!profilesfile.empty() and profilesfile != "Null") {
            if (flatfitness) {
                std::cerr << "Giving a preferences profiles file and setting 'flatfitness' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            prefs = open_preferences(profilesfile, Nstate);
            loaded_profiles = true;
            Ncat = static_cast<int>(std::get<0>(prefs).size());
            if (static_cast<int>(std::get<1>(prefs).size()) == Nsite) {
                std::cout << "Input preferences profiles file has the same number of profiles "
                             "as the alignment, hence site allocations are considered fixed."
                          << std::endl;
                loaded_profiles_allocation = true;

            } else {
                std::cout << "Input preferences profiles file has a smaller number of profiles "
                             "than the alignment, hence site allocations are computed."
                          << std::endl;
            }
        }
        // explicit checks (asserts vanish in Release): profiles read from file may exceed Nsite
        if (Ncat < 1 or Ncat > Nsite) {
            std::cerr << "error: number of profile components (" << Ncat
                      << ") must be between 1 and the number of sites (" << Nsite << ")\n";
            exit(1);
        }

        std::vector<double> opened_delta_omega_array;
        if (!delta_omega_array_file.empty() and delta_omega_array_file != "Null") {
            opened_delta_omega_array = open_delta_omega_array(delta_omega_array_file);
            loaded_delta_omega_array = true;
            std::cout << "Input omega values are provided, hence are considered fixed but site "
                         "allocations are computed."
                      << std::endl;
            omegamode = 0;
            omegaNcat = opened_delta_omega_array.size();
        }
        // validate before any allocation depending on omegaNcat
        if (omegaNcat < 1) {
            std::cerr << "error: omegancat must be >= 1\n";
            exit(1);
        }
        if (omega_shift < 0) {
            std::cerr << "error: omegashift must be >= 0\n";
            exit(1);
        }
        if (omegaNcat > 1 and omegamode == 3) {
            std::cerr << "warning: omegancat = " << omegaNcat
                      << " without --freeomega: all omega components are fixed to 1, so the "
                         "components are identical\n";
        }
        std::cout << "Number of sites: " << Nsite << std::endl;
        std::cout << "Number of profile categories (cut-off for infinite mixture): " << Ncat
                  << '\n';
        std::cout << "Number of omega categories (finite mixture): " << omegaNcat << '\n';
        std::cout << "Additive shift in omega: " << omega_shift << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        if (!tree_stream) {
            std::cerr << "error: cannot open tree file " << treefile << std::endl;
            exit(1);
        }
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
        nucrelratehyperinvconc = 0.1;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix

        std::vector<double> nucstat_{};
        if (!nucstatfile.empty() and nucstatfile != "Null") {
            if (flatnucstat) {
                std::cerr << "Giving a nucstat file and setting 'flatnucstat' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            nucstat_ = open_nucstat(nucstatfile, Nnuc);
            loaded_nucstat = true;
        }
        std::vector<double> nucrelrate_{};
        if (!nucrelratefile.empty() and nucrelratefile != "Null") {
            if (flatnucrelrate) {
                std::cerr << "Giving a nucrelrate file and setting 'flatnucrelrate' to true "
                             "is incompatible. These parameters are exclusive to one another."
                          << std::endl;
                exit(1);
            }
            nucrelrate_ = open_nucrelrate(nucrelratefile, Nrr);
            loaded_nucrelrate = true;
        }
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, std::vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        if (loaded_nucrelrate) {
            std::copy(nucrelrate_.begin(), nucrelrate_.end(), nucrelrate.begin());
            clamp_nucrelrate = true;
        } else if (flatnucrelrate) {
            nucrelrate.assign(Nrr, 1.0 / Nrr);
            clamp_nucrelrate = true;
        }
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, std::vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        if (loaded_nucstat) {
            std::copy(nucstat_.begin(), nucstat_.end(), nucstat.begin());
            clamp_nucstat = true;
        } else if (flatnucstat) {
            nucstat.assign(Nnuc, 1.0 / Nnuc);
            clamp_nucstat = true;
        }

        nucmatrix = new GTRSubMatrixCodon(Nnuc, nucrelrate, nucstat, true);

        // base distribution (can be skipped)
        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat, basekappa);
        baseprofile_occupancy = new OccupancySuffStat(baseNcat);

        basecenterhypercenter.assign(
            Nstate, 1.0 / Nstate);
        basecenterhyperinvconc = 1.0 / Nstate;

        basecenterarray = new IIDDirichlet(baseNcat, basecenterhypercenter, basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Nstate;
        baseconchyperinvshape = 1.0;
        baseconcentrationarray = new IIDGamma(baseNcat, baseconchypermean, baseconchyperinvshape);
        for (int k = 0; k < baseNcat; k++) {
            (*baseconcentrationarray)[k] = Nstate;
        }
        if (basemin == 1) { (*baseconcentrationarray)[0] = 1.0; }
        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat, Nstate);
        componentalloc = new MultinomialAllocationVector(Ncat, baseweight->GetArray());
        componentcenterarray =
            new MixtureSelector<std::vector<double>>(basecenterarray, componentalloc);
        componentconcentrationarray =
            new MixtureSelector<double>(baseconcentrationarray, componentalloc);

        //
        // (truncated) Dirichlet mixture of aa fitness profiles
        //

        // Ncat fitness profiles iid from the base distribution
        componentcodonfitnessarray =
            new MultiDirichlet(componentcenterarray, componentconcentrationarray);
        if (flatfitness) { 
            componentcodonfitnessarray->Flatten();  
            clamp_profiles = true;
            clamp_profiles_allocation = true;
        }

        if(Ncat == 1) {
            clamp_profiles_allocation = true;
        }

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
        if (omegamode == 3 and !loaded_delta_omega_array) {
            // omega fixed: pin all components to 1.0 (not just the first), so a
            // finite mixture with omegancat>1 does not run with frozen random draws.
            for (int omegacat = 0; omegacat < omegaNcat; omegacat++) {
                (*delta_omega_array)[omegacat] = 1.0;
            }
        }
        if (loaded_delta_omega_array) {
            for (int omegacat = 0; omegacat < omegaNcat; omegacat++) {
                (*delta_omega_array)[omegacat] = opened_delta_omega_array[omegacat];
            }
            clamp_delta_omega_array = true;
        }
        // will be a constant
        omega_weight_kappa = 1.0;
        omega_weight = new Dirichlet(omegaNcat, omega_weight_kappa);
        omega_alloc = new MultinomialAllocationVector(GetNsite(), omega_weight->GetArray());

        // selector, specifying which aa fitness array should be used for each site
        sitecodonfitnessarray =
            new MixtureSelector<std::vector<double>>(componentcodonfitnessarray, profile_alloc);

        if (loaded_profiles) {
            for (int cat = 0; cat < Ncat; cat++) {
                (*componentcodonfitnessarray)[cat] = std::get<0>(prefs)[cat];
            }
            clamp_profiles = true;
        }
        if (loaded_profiles_allocation) {
            for (int site = 0; site < Nsite; site++) {
                (*profile_alloc)[site] = static_cast<unsigned>(std::get<1>(prefs)[site]);
                assert(sitecodonfitnessarray->GetVal(site).size() ==
                       Nstate);
            }
            clamp_profiles_allocation = true;
        }

        // Ncat*omegaNcat mut sel codon matrices (based on the Ncat fitness profiles of the mixture,
        // and the omegaNcat finite mixture for omega)
        componentcodonmatrixbidimarray = new CodonMutSelCodonMatrixBidimArray(GetCodonStateSpace(),
            nucmatrix, componentcodonfitnessarray, delta_omega_array, omega_shift);

        // Bidimselector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(
            componentcodonmatrixbidimarray, profile_alloc, omega_alloc);

        // Bidimselector, specifying which codon matrix should be used for each site (needed to
        // collect omegapathsuffstatarray)
        sitecodonsubmatrixarray = new DoubleMixtureSelector<CodonMutSelOmegaCodonSubMatrix>(
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
        model_node(info, "componentcodonfitnessarray", *componentcodonfitnessarray);
        model_node(info, "profile_alloc", *profile_alloc);
        model_node(info, "delta_omega_array", *delta_omega_array);
        model_node(info, "omega_alloc", *omega_alloc);
        model_node(info, "omega_weight", *omega_weight);
        model_node(info, "delta_omegahyperinvshape", delta_omegahyperinvshape);
        model_node(info, "delta_omegahypermean", delta_omegahypermean);


        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        model_stat(info, "length", [this]() { return branchlength->GetTotalLength(); });
        model_stat(info, "omega", [this]() { return GetMeanComponentOmega(); });
        model_stat(info, "ds", [this]() { return GetPredictedRelativedS(); });
        model_stat(info, "dn", [this]() { return GetPredictedRelativedN(); });
        model_stat(info, "dnds", [this]() { return GetPredictedRelativedNdS(); });
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
        // kappa/basekappa are declared both as model nodes (.chain file) and as stats (.trace
        // file): these are two distinct files, so each column appears once per file.
        model_stat(info, "kappa", kappa);
        if (baseNcat > 1) {
            model_stat(info, "basencluster", [this]() { return GetBaseNcluster(); });
            if (basemin) {
                model_stat(info, "baseweight1", [this]() { return baseweight->GetVal(1); });
            } else {
                model_stat(info, "basekappa", basekappa);
            }
        }
        model_stat(info, "codonent", [this]() { return GetMeanEntropy(); });
        model_stat(info, "meancodonconc", [this]() { return GetMeanComponentConcentration(); });
        model_stat(info, "codoncenterent", [this]() { return GetMeanComponentCodonEntropy(); });
        model_stat(info, "rrent", [this]() { return GetNucRREntropy(); });
        model_stat(info, "statent", [this]() { return GetNucStatEntropy(); });
    }
    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    bool FlatFitness() const { return flatfitness; }
    bool FlatNucStat() const { return flatnucstat; }
    bool FlatNucRelRate() const { return flatnucrelrate; }

    //! return branch length
    double GetBranchLength(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return branchlength->GetVal(tree->branch_index(node));
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! return current omega value for omega mixture of component k
    double GetComponentOmega(int k) const { return omega_shift + delta_omega_array->GetVal(k); }

    double GetNucRR(int i) { return nucrelrate[i]; }
    double GetNucStat(int i) { return nucstat[i]; }
    double GetAC() {return nucrelrate[0];}
    double GetAG() {return nucrelrate[1];}
    double GetAT() {return nucrelrate[2];}
    double GetCG() {return nucrelrate[3];}
    double GetCT() {return nucrelrate[4];}
    double GetGT() {return nucrelrate[5];}
    double GetA() { return nucstat[0]; }
    double GetC() { return nucstat[1]; }
    double GetG() { return nucstat[2]; }
    double GetT() { return nucstat[3]; }


    double GetSiteCodonFitness(int site, int i) {
        // sitecodonfitnessarray is a MixtureSelector that already applies
        // profile_alloc in GetVal(site); do NOT index it again with
        // profile_alloc->GetVal(site) (that double-applies the allocation and
        // returns the wrong profile whenever Ncat>1 and the site is not in
        // profile 0).
        return sitecodonfitnessarray->GetVal(site)[i];
    }
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

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
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
        if (nucmode < 2) {
            if (!clamp_nucstat) { total += NucStatLogPrior(); }
            if (!clamp_nucrelrate) { total += NucRelRateLogPrior(); }
        }
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
            total += CodonLogPrior();
            // allocations given the mixture weights
            total += profile_alloc->GetLogProb();
            total += componentalloc->GetLogProb();
        }
        if (omegamode < 2) { total += DeltaOmegaLogProb(); }
        // omega mixture: Dirichlet prior on the weights + allocations given the weights
        total += omega_weight->GetLogProb();
        total += omega_alloc->GetLogProb();
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


    //! log prior over nuc rates rho and pi (uniform)
    double NucStatLogPrior() const {
        double total = 0;
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }
    double NucRelRateLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        return total;
    }
    //! log prior over concentration parameters basekappa of mixture of base
    //! distribution
    double BaseStickBreakingHyperLogPrior() const { return -basekappa / 10; }

    //! log prior over weights of stick breaking process of base distribution
    double BaseStickBreakingLogPrior() const { return baseweight->GetLogProb(basekappa); }

    //! log prior over concentration parameters kappa of stick-breaking mixture of
    //! codon fitness profiles
    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    //! log prior over weights of stick breaking process of codon fitness profiles
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

    //! log prior of codon fitness profiles
    double CodonLogPrior() const { return componentcodonfitnessarray->GetLogProb(); }

    //! log prior of codon fitness profile k
    double CodonLogPrior(int k) const { return componentcodonfitnessarray->GetLogProb(k); }

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
        for (int site : profile_sites[k]) {
            tot += sitepathsuffstatarray->GetVal(site).GetLogProb(sitesubmatrixarray->GetVal(site));
        }
        return tot;
    }

    //! return log prob of the substitution mappings over sites allocated to omega
    //! component k of the omega mixture
    double PathSuffStatOmegaLogProb(int k) const {
        return componentomegapathsuffstatarray->GetVal(k).GetLogProb(GetComponentOmega(k));
    }

    //! return log prob of first-level mixture components (i.e. all codon
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(
            basecenterarray->GetVal(k), 1.0 / baseconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------
    //! log prob factor to be recomputed when moving nucleotide mutation rate and codonfitness
    //! parameters (nucrelrate, nucstat and codonfitness)
    
    double NucRatesCodonLogProb() const { return  NucRelRateLogPrior() + NucStatLogPrior() + CodonLogPrior() + PathSuffStatLogProb(); }

    double NucRatesCodonLogProb(int k) const { return NucRelRateLogPrior() + NucStatLogPrior() + CodonLogPrior(k) + PathSuffStatLogProb(k); }

    //! log prob factor to be recomputed when moving nucleotide mutation rate
    //! parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

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
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete MCMC move schedule on all parameters
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            
            if (blmode < 2) { 
                MoveBranchLengths();
                //MoveBranchLengthsGlobal(0.1, 5);
            }
            
            CollectSitePathSuffStat();

            if (nucmode < 2) { 
                
                if (!clamp_nucrelrate) {
                    MoveNucRates();
                }
                if (!clamp_nucstat) {
                    MoveNucStat(); 
                }
            }
           
            if (!clamp_profiles || !clamp_profiles_allocation) {
                profilechrono.Start();
                MoveCodonMixture(3);
                profilechrono.Stop();

                basechrono.Start();
                if (basemode < 2 and !clamp_profiles) { MoveBase(3); }
                basechrono.Stop();
            }
               
            if (omegamode < 2) { MoveOmegaMixture(1); }

  
            totchrono.Stop();
        }
    }
    
    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1 && Ncat > 1) { ResampleBaseAlloc(); }
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

    //! acceptance rate of a move; 0 when no proposal was made (e.g. no occupied component)
    static double AcceptanceRate(double nacc, double ntot) { return ntot > 0 ? nacc / ntot : 0.0; }

    //! MH move on nucleotide stationary probabilities.
    //! Cost note: nucleotide parameters affect every codon matrix, so each proposal
    //! rebuilds all matrices (UpdateMatrices) and rescans all sites; this is inherent.
    void MoveNucStat() {
        for (double tuning : {0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01}) {
            Move::Profile(nucstat, tuning, 1, 1, &CodonMutSelMultipleOmegaModel::NucRatesLogProb,
                &CodonMutSelMultipleOmegaModel::UpdateMatrices, this);
        }
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        for (double tuning : {0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01}) {
            Move::Profile(nucrelrate, tuning, 2, 1, &CodonMutSelMultipleOmegaModel::NucRatesLogProb,
                &CodonMutSelMultipleOmegaModel::UpdateMatrices, this);
            Move::Profile(nucrelrate, tuning, 1, 1, &CodonMutSelMultipleOmegaModel::NucRatesLogProb,
                &CodonMutSelMultipleOmegaModel::UpdateMatrices, this);
        }
    }



    void MoveCodonMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {           
            if (!clamp_profiles) {
                MoveCodonProfiles();
                ResampleEmptyProfileComponents();
            }
            if (!clamp_profiles_allocation) {
                ResampleProfileAlloc();
                ProfileLabelSwitchingMove();
                ResampleProfileWeights();
                MoveKappa();
            }
            if (!clamp_profiles_allocation || !clamp_profiles) { CorruptCodonMatrices(); }
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyProfileComponents() {
        componentcodonfitnessarray->PriorResample(*profile_occupancy);
        componentcodonmatrixbidimarray->CorruptCodonMatricesRowOccupancy(*profile_occupancy);
    }

    //! MH move on codon fitness profiles (occupied components only)
    void MoveCodonProfiles() {
        CompMoveCodonProfiles(1);
        MulMoveCodonProfiles(1);
    }


    //! MH move on codon fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveCodonProfiles(int nrep) {
        MoveCodonGamma(3.0, nrep);
        MoveCodonGamma(1.0, nrep);
        return 1.0;
    }

    //! MH move on codon fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveCodon(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(Nstate);
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                std::vector<double> &codon = (*componentcodonfitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Nstate; l++) {
                        bk[l] = codon[l];
                    }
                    double deltalogprob = -CodonLogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(
                        codon, Nstate, tuning, n);
                    deltalogprob += loghastings;
                    CorruptProfileCodonMatrices(i);
                    deltalogprob += CodonLogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Nstate; l++) {
                            codon[l] = bk[l];
                        }
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return AcceptanceRate(nacc, ntot);
    }


    //! MH move on codon fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveCodonProfiles(int nrep) {
        MoveCodon(0.3, 2, 1);
        MoveCodon(0.15, 4, 1);
        MoveCodon(0.07, 8, 1);
        MoveCodon(0.03, 16, 1);
        MoveCodon(0.01, 32, 1);
        MoveCodon(1, 2, 1);
        
        return 1.0;
    }
        
    
    //! helper function: log density of Nstate gammas
    double GammaCodonLogPrior(
        const std::vector<double> &x, const std::vector<double> &codoncenter, double codonconc) {
        double total = 0;
        for (int l = 0; l < Nstate; l++) {
            total += (codonconc * codoncenter[l] - 1) * log(x[l]) - x[l] -
                     Random::logGamma(codonconc * codoncenter[l]);
        }
        return total;
    }

    //! MH move on codon fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MoveCodonGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Ncat; i++) {
            if (profile_occupancy->GetVal(i)) {
                double codonconc = componentconcentrationarray->GetVal(i);
                const std::vector<double> &codoncenter = componentcenterarray->GetVal(i);

                std::vector<double> &codon = (*componentcodonfitnessarray)[i];
                std::vector<double> x(Nstate, 0);
                double z = Random::sGamma(codonconc);
                for (int l = 0; l < Nstate; l++) { x[l] = z * codon[l]; }

                double bkz = z;
                std::vector<double> bkx = x;
                std::vector<double> bkcodon = codon;

                for (int rep = 0; rep < nrep; rep++) {
                    double deltalogprob =
                        -GammaCodonLogPrior(x, codoncenter, codonconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l = 0; l < Nstate; l++) {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l = 0; l < Nstate; l++) {
                        codon[l] = x[l] / z;
                        if (codon[l] < 1e-50) { codon[l] = 1e-50; }
                    }

                    deltalogprob += loghastings;

                    CorruptProfileCodonMatrices(i);

                    deltalogprob +=
                        GammaCodonLogPrior(x, codoncenter, codonconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        bkcodon = codon;
                        bkx = x;
                        bkz = z;
                    } else {
                        codon = bkcodon;
                        x = bkx;
                        z = bkz;
                        CorruptProfileCodonMatrices(i);
                    }
                    ntot++;
                }
            }
        }
        return AcceptanceRate(nacc, ntot);
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
        UpdateProfileSites();
    }

    //! rebuild the list of sites allocated to each profile component
    void UpdateProfileSites() {
        profile_sites.assign(Ncat, std::vector<int>());
        for (int site = 0; site < Nsite; site++) {
            profile_sites[profile_alloc->GetVal(site)].push_back(site);
        }
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
        componentcodonfitnessarray->Permute(permut);
        UpdateProfileSites();
    }

    //! Gibbs resample mixture weights (based on profile_occupancy suff stats)
    void ResampleProfileWeights() { weight->GibbsResample(*profile_occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        Move::Scaling(kappa, 1.0, 10, &CodonMutSelMultipleOmegaModel::StickBreakingHyperLogProb,
            &CodonMutSelMultipleOmegaModel::NoUpdate, this);
        Move::Scaling(kappa, 0.3, 10, &CodonMutSelMultipleOmegaModel::StickBreakingHyperLogProb,
            &CodonMutSelMultipleOmegaModel::NoUpdate, this);
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
        componentcodonfitnessarray->AddSuffStat(*basesuffstatarray, *componentalloc);
    }

    //! MCMC module for moving the center parameters of the components of the the
    //! base mixture
    double MoveBaseCenters(double tuning, int n) {
        double nacc = 0;
        double ntot = 0;
        std::vector<double> bk(Nstate, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseprofile_occupancy->GetVal(k)) {
                std::vector<double> &codon = (*basecenterarray)[k];
                bk = codon;
                double deltalogprob = -BaseLogProb(k);
                double loghastings =
                    Random::ProfileProposeMove(codon, Nstate, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    codon = bk;
                }
                ntot++;
            }
        }
        return AcceptanceRate(nacc, ntot);
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
        return AcceptanceRate(nacc, ntot);
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
            double tmp = Random::logDirichletDensity(componentcodonfitnessarray->GetVal(cat),
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
            &CodonMutSelMultipleOmegaModel::BaseStickBreakingHyperLogProb,
            &CodonMutSelMultipleOmegaModel::NoUpdate, this);
        Move::Scaling(basekappa, 0.3, 10,
            &CodonMutSelMultipleOmegaModel::BaseStickBreakingHyperLogProb,
            &CodonMutSelMultipleOmegaModel::NoUpdate, this);
    }

    //! MCMC module for the mixture of omega
    void MoveOmegaMixture(int nrep) {
        CollectSiteOmegaPathSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            if (omegaNcat > 1) { 
                ResampleOmegaAlloc(); 
                UpdateOmegaOccupancies();
            }
            if (!clamp_delta_omega_array) {
                CollectComponentOmegaPathSuffStat();
                MoveOmegaValues(1.0, 3);
                MoveOmegaValues(0.1, 3);
                ResampleEmptyOmegaComponents();
                MoveOmegaHyper();
            }
            if (omegaNcat > 1) {
                ResampleOmegaWeights();
            }
        }
        // only omega changed here (profiles did not): corrupt per column, which
        // skips the fitness recompute (CorruptColCodonMatrices uses
        // SetOmegaNoCorrupt + CorruptMatrixNoFitnessRecomput).
        for (int j = 0; j < omegaNcat; j++) { componentcodonmatrixbidimarray->CorruptColCodonMatrices(j); }
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
        return AcceptanceRate(nacc, ntot);
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
                &CodonMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &CodonMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahypermean, 0.3, 10,
                &CodonMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &CodonMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 1.0, 10,
                &CodonMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &CodonMutSelMultipleOmegaModel::NoUpdate, this);
            Move::Scaling(delta_omegahyperinvshape, 0.3, 10,
                &CodonMutSelMultipleOmegaModel::DeltaOmegaLogProb,
                &CodonMutSelMultipleOmegaModel::NoUpdate, this);
        }
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return mutation rate between nucleotides from the mutation matrix
    double GetNucRate(int i, int j) const {
        nucmatrix->UpdateMatrix();
        return nucmatrix->RelativeRate(i, j);
    }


    //! return number of occupied components in first-level mixture (mixture of
    //! codon fitness profiles)
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


    //! return mean omega, weighted by the mixture weights (so empty components
    //! resampled from the prior do not bias the traced value)
    double GetMeanComponentOmega() const {
        double tot = 0;
        double totw = 0;
        for (int i{0}; i < omegaNcat; i++) {
            tot += omega_weight->GetVal(i) * GetComponentOmega(i);
            totw += omega_weight->GetVal(i);
        }
        return (totw > 0) ? tot / totw : 0.0;
    }

    //! return mean entropy of codon fitness profiles
    double GetMeanEntropy() const { return componentcodonfitnessarray->GetMeanEntropy(); }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentConcentration() const {
        double tot = 0;
        double totw = 0;
        for (int i = basemin; i < baseNcat; i++) {
            tot += baseprofile_occupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
            totw += baseprofile_occupancy->GetVal(i);
        }
        return (totw > 0) ? tot / totw : 0.0;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentCodonEntropy() const {
        double tot = 0;
        double totw = 0;
        for (int i = 0; i < baseNcat; i++) {
            tot +=
                baseprofile_occupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
            totw += baseprofile_occupancy->GetVal(i);
        }
        return (totw > 0) ? tot / totw : 0.0;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucstat); }

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

    double GetPredictedEffectivedN() const {
        double mean = 0;
        for (int i = 0; i < GetNsite(); i++) {
            mean += GetSiteOmega(i) * sitecodonsubmatrixarray->GetVal(i).GetPredictedDN();
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

    //! compute the relative dN, dS and dN/dS in a single pass over sites, so the
    //! three trace stats do not each recompute the full per-site flow.
    std::tuple<double, double, double> GetPredictedRelativeFlows() const {
        double totdn = 0, totds = 0;
        for (int i = 0; i < GetNsite(); i++) {
            double dn = 0, ds = 0;
            std::tie(dn, ds) = sitecodonsubmatrixarray->GetVal(i).GetRelativeFlowDNDS();
            // preserve original semantics: dN and dN/dS are scaled by site omega,
            // dS is not.
            totdn += GetSiteOmega(i) * dn;
            totds += ds;
        }
        double dn = totdn / GetNsite();
        double ds = totds / GetNsite();
        double dnds = (ds > 0) ? dn / ds : 0.0;
        return std::make_tuple(dn, ds, dnds);
    }

    double GetPredictedRelativedNdS() const {
        double dn = 0, ds = 0, dnds = 0;
        std::tie(dn, ds, dnds) = GetPredictedRelativeFlows();
        return dnds;
    }

    double GetPredictedRelativedN() const {
        double dn = 0, ds = 0, dnds = 0;
        std::tie(dn, ds, dnds) = GetPredictedRelativeFlows();
        return dn;
    }

    double GetPredictedRelativedS() const {
        double dn = 0, ds = 0, dnds = 0;
        std::tie(dn, ds, dnds) = GetPredictedRelativeFlows();
        return ds;
    }
    const Tree &GetTree() const { return *tree; }
    const int GetNcat() const { return Ncat; }
    const int GetOmegaNcat() const { return omegaNcat; }
    const std::vector<double> &GetProfile(int i) const { return sitecodonfitnessarray->GetVal(i); }
    const double GetProfileCodon(int i, int c) const {
        return componentcodonfitnessarray->GetVal(i)[c];
    }
    const int GetProfileAlloc(int i) const { return profile_alloc->GetVal(i); }
    const double GetOmega(int i) const { return delta_omega_array->GetVal(i); }
    const int GetOmegaAlloc(int i) const { return omega_alloc->GetVal(i); }
    static std::string GetModelName() { return CLASS_NAME(CodonMutSelMultipleOmega); }

    void GetModelStamp(std::ostream &os) const {
        os << GetModelName() << '\t';
        os << datafile << '\t' << treefile << '\t' << profilesfile << '\t';
        os << delta_omega_array_file << '\t';
        os << omegaNcat << '\t';
        os << omega_shift << '\t';
        os << (basemin ? -baseNcat : baseNcat) << '\t';
        os << Ncat << '\t';
        os << omegamode << '\t';
        os << flatfitness << '\t';
        os << flatnucstat << '\t';
        os << flatnucrelrate << '\t';
        // nucstatfile/nucrelratefile are appended at the END so the first 13
        // tokens (read by LikelihoodFreePhylogenetics' readBayescodeParametersMutSelC)
        // keep their original positions.
        os << nucstatfile << '\t' << nucrelratefile << '\t';
    }

    //! stamp with exactly 13 tokens, the format read by
    //! LikelihoodFreePhylogenetics' readBayescodeParametersMutSelC (Ncat at
    //! index 8). Used for --for_simulation .pvalues output. The full stamp
    //! (with nucstatfile/nucrelratefile appended) is written to .param for resume.
    void GetPvaluesStamp(std::ostream &os) const {
        os << GetModelName() << '\t';
        os << datafile << '\t' << treefile << '\t' << profilesfile << '\t';
        os << delta_omega_array_file << '\t';
        os << omegaNcat << '\t';
        os << omega_shift << '\t';
        os << (basemin ? -baseNcat : baseNcat) << '\t';
        os << Ncat << '\t';
        os << omegamode << '\t';
        os << flatfitness << '\t';
        os << flatnucstat << '\t';
        os << flatnucrelrate << '\t';
    }

    void ToStream(std::ostream &os) const {
        GetModelStamp(os);
        tracer->write_line(os);
    }

    CodonMutSelMultipleOmegaModel(std::istream &is) {
        std::string model_name;
        is >> model_name;
        if (model_name != GetModelName()) {
            std::cerr << "Expected " << GetModelName() << " for model name, got " << model_name
                      << "\n";
            exit(1);
        }
        is >> datafile >> treefile >> profilesfile;
        is >> delta_omega_array_file;
        is >> omegaNcat;
        is >> omega_shift;
        is >> baseNcat;
        is >> Ncat;
        is >> omegamode;
        is >> flatfitness;
        is >> flatnucstat;
        is >> flatnucrelrate;
        is >> nucstatfile >> nucrelratefile;
        init();
        tracer->read_line(is);
        Update();
    };
};