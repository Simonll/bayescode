
#include "CodonM8Model.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"

#include "Chrono.hpp"


class MultiGeneCodonM8Model : public MultiGeneMPIModule	{

    private:

	Tree* tree;
	FileSequenceAlignment* refdata;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    double purifmeanhypermean;
	double purifmeanhyperinvconc;
	IIDBeta* purifmeanarray;
	BetaSuffStat purifmeansuffstat;

    double purifinvconchypermean;
    double purifinvconchyperinvshape;
	IIDGamma* purifinvconcarray;
	GammaSuffStat purifinvconcsuffstat;

    double dposomhypermean;
    double dposomhyperinvshape;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    double poswhypermean;
    double poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    IIDDirichlet* purifweightarray;
    DirichletSuffStat purifweightsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double pi;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    NucPathSuffStat nucpathsuffstat;

    std::vector<CodonM8Model*> geneprocess;

    double* lnL;

    int ncat;

    int burnin;

    Chrono timepercycle;

    public:

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    MultiGeneCodonM8Model(string datafile, string intreefile, int inncat, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs), purifweightsuffstat(3) {

        burnin = 0;

        ncat = inncat;
        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;

        AllocateAlignments(datafile);
        treefile = intreefile;

        // all datafiles have all taxa (with missing data if needed) in same order
        // makes it easier to register tree with data, etc.

        string filename = genename[0];
        refdata = new FileSequenceAlignment(filename);
		refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        std::cerr << "number of taxa : " << Ntaxa << '\n';
        std::cerr << "number of branches : " << Nbranch << '\n';
        std::cerr << "-- Tree and data fit together\n";

        // Allocate();
    }

    void Allocate() {

        
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);
        lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelrate.assign(Nrr,0);
        double totrr = 0;
        for (int k=0; k<Nrr; k++)	{
            nucrelrate[k] = Random::sExpo();
            totrr += nucrelrate[k];
        }
        for (int k=0; k<Nrr; k++)	{
            nucrelrate[k] /= totrr;
        }

        nucstat.assign(Nnuc,0);
        double totstat = 0;
        for (int k=0; k<Nnuc; k++)	{
            nucstat[k] = Random::sGamma(1.0);
            totstat += nucstat[k];
        }
        for (int k=0; k<Nnuc; k++)	{
            nucstat[k] /= totstat;
        }
        nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        purifmeanhypermean = 0.5;
        purifmeanhyperinvconc = 0.5;
        purifinvconchypermean = 2.0;
        purifinvconchyperinvshape = 1.0;

        dposomhypermean = 1.0;
        dposomhyperinvshape = 1.0;

        pi = pihypermean;
        if (! pi)   {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        else    {
            poswhypermean = 0.1;
            poswhyperinvconc = 1;
        }

        purifweighthypercenter.assign(3,1.0/3);
        purifweighthyperinvconc = 1.0/3;

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        purifmeanarray = new IIDBeta(GetNgene(),purifmeanalpha,purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        purifinvconcarray = new IIDGamma(GetNgene(),purifinvconcalpha,purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
		dposomarray = new IIDGamma(GetNgene(),dposomalpha,dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetNgene(),pi,poswalpha,poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        purifweightarray = new IIDDirichlet(GetNgene(),purifweighthypercenter,purifweighthyperconc);

        lnL = new double[Ngene];

        cerr << "gene processes\n";
        if (! GetMyid())    {
            geneprocess.assign(0,(CodonM8Model*) 0);
        }
        else    {
            geneprocess.assign(Ngene,(CodonM8Model*) 0);

            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene] == myid)    {
                    geneprocess[gene] = new CodonM8Model(genename[gene],treefile,ncat,pi);
                }
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendHyperParameters();
            MasterSendGlobalParameters();
            MasterSendMixture();
        }
        else    {
            SlaveReceiveHyperParameters();
            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene] == myid)    {
                    geneprocess[gene]->Allocate();
                }
            }

            SlaveReceiveGlobalParameters();
            SlaveReceiveMixture();

            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene] == myid)    {
                    geneprocess[gene]->UpdateMatrices();
                    geneprocess[gene]->Unfold();
                }
            }
        }
    }

    void TraceHeader(ostream& os)   {

        os << "#time\t";
        os << "logprior\tlnL\tlength\t";
        os << "pi\t";
        os << "nposfrac\t";
        os << "meanposfrac\t";
        os << "meanposom\t";
        os << "meanw0\tmeanw1\tpurifinvconc\t";
        os << "purifmean\tinvconc\tpurifinvconc\tinvshape\t";
        os << "poswmean\tinvconc\t";
        os << "dposommean\tinvshape\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void MasterTrace(ostream& os)    {
        os << timepercycle.GetTime() << '\t';
        timepercycle.Stop();
        timepercycle.Reset();
        timepercycle.Start();
		os << GetLogPrior() << '\t';
        MasterReceiveLogLikelihood();
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
        os << pi << '\t';
        os << GetNpos() << '\t';
        os << GetMeanPosFrac() << '\t';
        os << GetMeanPosOmega() << '\t';
        os << purifweighthypercenter[0] << '\t' << purifweighthypercenter[2] << '\t' << purifweighthyperinvconc << '\t';
        os << purifmeanhypermean << '\t' << purifmeanhyperinvconc << '\t' << purifinvconchypermean << '\t' << purifinvconchyperinvshape << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
		os.flush();
    }

    void TracePosWeight(ostream& os) {

        for (int gene=0; gene<Ngene; gene++)    {
            os << (*poswarray)[gene] << '\t';
        }
        os << '\n';
        os.flush();
    }

    int GetNpos()    {
        return GetNgene() - poswarray->GetNullSet();
    }

    double GetMeanPosFrac() {
        return poswarray->GetPosMean();
    }

    double GetMeanPosOmega()    {

        double tot = 0;
        double totweight = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (poswarray->GetVal(gene)) {
                tot += poswarray->GetVal(gene) * dposomarray->GetVal(gene);
                totweight += poswarray->GetVal(gene);
            }
        }
        return tot / totweight + 1;
    }

    void SlaveTrace()   {
        SlaveSendLogLikelihood();
    }

    double GetLogPrior()    {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += HyperLogPrior();
        total += NucLogPrior();
		return total;
    }

    double NucLogPrior()    {
        return 0;
    }

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double HyperLogPrior()   {
        double total = 0;
        total -= purifmeanhyperinvconc;
        total -= purifinvconchypermean;
        total -= purifinvconchyperinvshape;
        total -= poswhyperinvconc;
        total -= dposomhypermean;
        total -= dposomhyperinvshape;
        total -= purifweighthyperinvconc;
        return total;
    }

    double HyperSuffStatLogProb()   {
        double total = 0;

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        total += purifmeansuffstat.GetLogProb(purifmeanalpha,purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        total += purifinvconcsuffstat.GetLogProb(purifinvconcalpha,purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha,dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pi,poswalpha,poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        total += purifweightsuffstat.GetLogProb(purifweighthypercenter,purifweighthyperconc);

        return total;
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    double GetLogLikelihood()   {
        double tot = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            tot += lnL[gene];
        }
        return tot;
    }

	double GetTotalLength()	{
		double tot = 0;
		for (int j=1; j<Nbranch; j++)	{
			tot += branchlength->GetVal(j);
		}
		return tot;
	}

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

	void Monitor(ostream& os) {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

    void MasterMove() {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveLengthSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalParameters();

            MasterReceiveMixture();
            MasterMoveMixtureHyperParameters();
            MasterSendHyperParameters();

            MasterReceiveNucPathSuffStat();
            MasterMoveNuc();
            MasterSendGlobalParameters();
        }
    }

    // slave move
    void SlaveMove() {

        SlaveResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalParameters();

            SlaveCollectPathSuffStat();
            SlaveMoveOmega();
            SlaveSendMixture();
            SlaveReceiveHyperParameters();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalParameters();

        }
    }

    void MasterSendGlobalParameters() {

        int N = Nbranch + Nrr + Nnuc;
        double* array = new double[N];
        int i = 0;
        for (int j=0; j<Nbranch; j++)   {
            array[i++] = branchlength->GetVal(j);
        }
        for (int j=0; j<Nrr; j++)   {
            array[i++] = nucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            array[i++] = nucstat[j];
        }
        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalParameters()   {

        int N = Nbranch + Nrr + Nnuc;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = array[i++];
        }
        for (int j=0; j<Nrr; j++)   {
            nucrelrate[j] = array[i++];
        }
        for (int j=0; j<Nnuc; j++)  {
            nucstat[j] = array[i++];
        }

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->SetBranchLengths(*branchlength);
                geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
            }
        }
        delete[] array;
    }

    void MasterSendHyperParameters() {

        int N = 13;
        double* array = new double[N];
        int i = 0;
        array[i++] = purifmeanhypermean;
        array[i++] = purifmeanhyperinvconc;
        array[i++] = purifinvconchypermean;
        array[i++] = purifinvconchyperinvshape;
        array[i++] = pi;
        array[i++] = poswhypermean;
        array[i++] = poswhyperinvconc;
        array[i++] = dposomhypermean;
        array[i++] = dposomhyperinvshape;
        array[i++] = purifweighthypercenter[0];
        array[i++] = purifweighthypercenter[1];
        array[i++] = purifweighthypercenter[2];
        array[i++] = purifweighthyperinvconc;

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveHyperParameters()   {

        int N = 13;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        purifmeanhypermean = array[i++];
        purifmeanhyperinvconc = array[i++];
        purifinvconchypermean = array[i++];
        purifinvconchyperinvshape = array[i++];
        pi = array[i++];
        poswhypermean = array[i++];
        poswhyperinvconc = array[i++];
        dposomhypermean = array[i++];
        dposomhyperinvshape = array[i++];
        purifweighthypercenter[0] = array[i++];
        purifweighthypercenter[1] = array[i++];
        purifweighthypercenter[2] = array[i++];
        purifweighthyperinvconc = array[i++];

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->SetMixtureHyperParameters(purifmeanhypermean,purifmeanhyperinvconc,purifinvconchypermean,purifinvconchyperinvshape,dposomhypermean,dposomhyperinvshape,pi,poswhypermean,poswhyperinvconc,purifweighthypercenter,purifweighthyperinvconc);
            }
        }
        delete[] array;
    }

    void SlaveResampleSub()  {

        double frac = 1.0;
        if (burnin > 10)    {
            frac = 0.2;
        }
        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->ResampleSub(frac);
            }
        }
        burnin++;
    }

    void SlaveSendLengthSuffStat()  {

        lengthsuffstatarray->Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectLengthSuffStat();
                lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
            }
        }

        int* count = new int[Nbranch];
        double* beta = new double[Nbranch];
        lengthsuffstatarray->Push(count,beta);
        MPI_Send(count,Nbranch,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Nbranch,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] count;
        delete[] beta;
    }

    void MasterReceiveLengthSuffStat()  {

        int* count = new int[Nbranch];
        double* beta = new double[Nbranch];
        lengthsuffstatarray->Clear();
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nbranch,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Nbranch,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            lengthsuffstatarray->Add(count,beta);
        }
        delete[] count;
        delete[] beta;
    }

    void MasterResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthsuffstatarray);
    }

	void MasterMoveLambda()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
		MoveLambda(1.0,10);
		MoveLambda(0.3,10);
		branchlength->SetScale(lambda);
    }

	double MoveLambda(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - LambdaLogProb() - LengthSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			lambda *= e;
			deltalogprob += LambdaLogProb() + LengthSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				lambda /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void MasterMoveMixtureHyperParameters()  {

        CollectHyperSuffStat();

		HyperSlidingMove(purifmeanhypermean,1.0,10,0,1);
		HyperSlidingMove(purifmeanhypermean,0.3,10,0,1);
		HyperScalingMove(purifmeanhyperinvconc,1.0,10);
		HyperScalingMove(purifmeanhyperinvconc,0.3,10);

		HyperScalingMove(purifinvconchypermean,1.0,10);
		HyperScalingMove(purifinvconchypermean,0.3,10);
		HyperScalingMove(purifinvconchyperinvshape,1.0,10);
		HyperScalingMove(purifinvconchyperinvshape,0.3,10);

		HyperSlidingMove(poswhypermean,1.0,10,0,1);
		HyperSlidingMove(poswhypermean,0.3,10,0,1);
		HyperScalingMove(poswhyperinvconc,1.0,10);
		HyperScalingMove(poswhyperinvconc,0.3,10);

		HyperScalingMove(dposomhypermean,1.0,10);
		HyperScalingMove(dposomhypermean,0.3,10);
		HyperScalingMove(dposomhyperinvshape,1.0,10);
		HyperScalingMove(dposomhyperinvshape,0.3,10);

        HyperProfileMove(purifweighthypercenter,1.0,1,10);
        HyperProfileMove(purifweighthypercenter,0.3,1,10);
        HyperScalingMove(purifweighthyperinvconc,1.0,10);
        HyperScalingMove(purifweighthyperinvconc,0.3,10);

        if (burnin > 10)    {
            if (pihyperinvconc)    {
                ResamplePi();
            }
        }

        SetArrays();
    }

    void CollectHyperSuffStat() {

		purifmeansuffstat.Clear();
		purifmeanarray->AddSuffStat(purifmeansuffstat);
		purifinvconcsuffstat.Clear();
		purifinvconcarray->AddSuffStat(purifinvconcsuffstat);
		poswsuffstat.Clear();
		poswarray->AddSuffStat(poswsuffstat);
		dposomsuffstat.Clear();
		dposomarray->AddSuffStat(dposomsuffstat);
        purifweightsuffstat.Clear();
        purifweightarray->AddSuffStat(purifweightsuffstat);
    }

    void SetArrays()    {

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        purifmeanarray->SetAlpha(purifmeanalpha);
        purifmeanarray->SetBeta(purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        purifinvconcarray->SetShape(purifinvconcalpha);
        purifinvconcarray->SetScale(purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        purifweightarray->SetCenter(purifweighthypercenter);
        purifweightarray->SetConcentration(purifweighthyperconc);
    }

    void ResamplePi()   {

        int n0 = poswsuffstat.GetN0();
        int n1 = poswsuffstat.GetN1();
        if ((n0+n1) != Ngene)   {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double pialpha = pihypermean / pihyperinvconc;
        double pibeta = (1-pihypermean) / pihyperinvconc;
        double a0 = Random::sGamma(pialpha + n0);
        double a1 = Random::sGamma(pibeta + n1);
        pi = a1 / (a0 + a1);
    }

	double HyperSlidingMove(double& x, double tuning, int nrep, double min = 0, double max = 0)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
		    x += m;
            if (max > min)  {
                while ((x < min) || (x > max))  {
                    if (x < min)    {
                        x = 2*min - x;
                    }
                    if (x > max)    {
                        x = 2*max - x;
                    }
                }
            }
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x -= m;
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double HyperScalingMove(double& x, double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
		    x *= e;
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double HyperProfileMove(vector<double>& x, double tuning, int n, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        vector<double> bk(x.size(),0);
		for (int rep=0; rep<nrep; rep++)	{
            bk = x;
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			deltalogprob += loghastings;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x = bk;
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void SlaveCollectPathSuffStat() {
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectPathSuffStat();
            }
        }
    }

    void SlaveSendNucPathSuffStat()  {

        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectComponentPathSuffStat();
                geneprocess[gene]->CollectNucPathSuffStat();
                nucpathsuffstat.Add(geneprocess[gene]->GetNucPathSuffStat());
            }
        }

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat.Push(count,beta);
        MPI_Send(count,Nnuc+Nnuc*Nnuc,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Nnuc*Nnuc,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] count;
        delete[] beta;
    }

    void MasterReceiveNucPathSuffStat()  {

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat.Clear();
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nnuc+Nnuc*Nnuc,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Nnuc*Nnuc,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            nucpathsuffstat.Add(count,beta);
        }
        delete[] count;
        delete[] beta;
    }

    void MasterMoveNuc()    {

		MoveRR(0.1,1,3);
		MoveRR(0.03,3,3);
		MoveRR(0.01,3,3);

		MoveNucStat(0.1,1,3);
		MoveNucStat(0.01,1,3);
    }

    double NucPathSuffStatLogProb() {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	double MoveRR(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nrr];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nrr; l++)	{
				bk[l] = nucrelrate[l];
			}
			double deltalogprob = -NucPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
			deltalogprob += loghastings;
            UpdateNucMatrix();
			deltalogprob += NucPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nrr; l++)	{
					nucrelrate[l] = bk[l];
				}
                UpdateNucMatrix();
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveNucStat(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nnuc];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nnuc; l++)	{
				bk[l] = nucstat[l];
			}
			double deltalogprob = -NucPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
			deltalogprob += loghastings;
            UpdateNucMatrix();
			deltalogprob += NucPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nnuc; l++)	{
					nucstat[l] = bk[l];
				}
                UpdateNucMatrix();
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void SlaveMoveOmega()  {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->MoveOmega();
            }
        }
    }

    void SlaveSendMixture()   {

        double* array = new double[7*Ngene];
        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                array[gene] = geneprocess[gene]->GetPurifMean();
                array[Ngene+gene] = geneprocess[gene]->GetPurifInvConc();
                array[2*Ngene+gene] = geneprocess[gene]->GetPosW();
                array[3*Ngene+gene] = geneprocess[gene]->GetDPosOm();
                array[4*Ngene+gene] = geneprocess[gene]->GetPurifWeight(0);
                array[5*Ngene+gene] = geneprocess[gene]->GetPurifWeight(1);
                array[6*Ngene+gene] = geneprocess[gene]->GetPurifWeight(2);
            }
            else    {
                array[gene] = -1;
                array[Ngene+gene] = -1;
                array[2*Ngene+gene] = -1;
                array[3*Ngene+gene] = -1;
                array[4*Ngene+gene] = -1;
                array[5*Ngene+gene] = -1;
                array[6*Ngene+gene] = -1;
            }
        }
        MPI_Send(array,7*Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveMixture()    {

        double* array = new double[7*Ngene];
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(array,7*Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<Ngene; gene++)    {
                if (array[gene] != -1)    {
                    (*purifmeanarray)[gene] = array[gene];
                    (*purifinvconcarray)[gene] = array[Ngene+gene];
                    (*poswarray)[gene] = array[2*Ngene+gene];
                    (*dposomarray)[gene] = array[3*Ngene+gene];
                    (*purifweightarray)[gene][0] = array[4*Ngene+gene];
                    (*purifweightarray)[gene][1] = array[5*Ngene+gene];
                    (*purifweightarray)[gene][2] = array[6*Ngene+gene];
                }
            }
        }
        delete[] array;
    }

    void MasterSendMixture()  {

        double* array = new double[7*Ngene];
        for (int gene=0; gene<Ngene; gene++)    {
            array[gene] = (*purifmeanarray)[gene];
            array[Ngene+gene] = (*purifinvconcarray)[gene];
            array[2*Ngene+gene] = (*poswarray)[gene];
            array[3*Ngene+gene] = (*dposomarray)[gene];
            array[4*Ngene+gene] = (*purifweightarray)[gene][0];
            array[5*Ngene+gene] = (*purifweightarray)[gene][1];
            array[6*Ngene+gene] = (*purifweightarray)[gene][2];
        }
        MPI_Bcast(array,7*Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveMixture()    {

        double* array = new double[7*Ngene];
        MPI_Bcast(array,7*Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);

        for (int gene=0; gene<Ngene; gene++)    {
            (*purifmeanarray)[gene] = array[gene];
            (*purifinvconcarray)[gene] = array[Ngene+gene];
            (*poswarray)[gene] = array[2*Ngene+gene];
            (*dposomarray)[gene] = array[3*Ngene+gene];
            (*purifweightarray)[gene][0] = array[4*Ngene+gene];
            (*purifweightarray)[gene][1] = array[5*Ngene+gene];
            (*purifweightarray)[gene][2] = array[6*Ngene+gene];
        }

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->SetMixtureParameters((*purifmeanarray)[gene],(*purifinvconcarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
            }
        }
        delete[] array;
    }

    void SlaveSendLogLikelihood()   {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                lnL[gene] = geneprocess[gene]->GetLogLikelihood();
            }
            else    {
                lnL[gene] = 0;
            }
        }
        MPI_Send(lnL,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    void MasterReceiveLogLikelihood()    {

        double* array = new double[Ngene];
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(array,Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<Ngene; gene++)    {
                if (array[gene])    {
                    lnL[gene] = array[gene];
                }
            }
        }
        delete[] array;
    }

};
