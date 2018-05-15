#include "mutation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <climits>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

extern MTRand rnd; // random number generator
extern FILE * fichierS; // output file


/*----------------------------------------------------------
Function recursion: iterates the life cycle.
Parameters are:
Nv: population size (nb of haploid organisms)
sigv: std deviation of distribution of mutational effects
nbSv: number of selected loci (affecting the phenotype) per genome; equal to number of cistronic regulators
NbGenv: number of generations
NbPrelimv: number of preliminary generations with recombination between sex chromosomes
pasv: time interval between measurements
s: mutational effect (mean) on allele
I: intensity of stabilizing selection
U_g: mutation rate for genes
U_c: mutation rate for cis regulators
Ut_dro: mutation rate for trans regulator: drosophila model
Ut_cel: mutation rate for trans regulator: C elegans model
Ut_mam: mutation rate for trans regulator: mammal model
Rg: genetic distance between genes
Rc: genetic distance between a cis-regulator and its regulated gene !!! RECOMBINATION FUNCTION ONLY VALID FOR Rc < Rg !!!!
 output: 0 if all loci recorded
-----------------------------------------------------------*/



void recursion(int Nv, double sigv, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double Ut_dro, double Ut_cel, double Ut_mam, double Rg, double Rc, int Rep, int output, double** allAverages)
{

	int i, j, k, nm, gen, mut, p1, p2, indiv, site, chrom, off_sex, Nmales_1, Nfemales_1, NjuvM, NjuvF;
	double w, wmmax, wfmax, h, dg, dc, dt, varm;

	double h0 = 0.25;
    double Q0 = 2;
    int N_1 = Nv - 1;
    int nbSm1 = nbSv - 1;
    int twonbSm1 = 2 * nbSv - 1;
    int Nmales = Nv / 2;
    double MLength = Rg * (nbSv - 1) + Rc;
    double UgTot = 2*Nv*U_g*nbSv;
    double UcTot = 2*Nv*U_c*nbSv;
    double UTot_dro = 2*Nv*Ut_dro;
    double UTot_cel = 2*Nv*Ut_cel;
    double UTot_mam = 2*Nv*Ut_mam;
    double expdom = -log(h0) / log(2.0);

    // creates result file (with parameter values in file name):
    char nomFichier[256];
    stringstream nomF;
    nomF << "N" << Nv << "_nbS" << nbSv << "_NbPrelim" << NbPrelimv << "_s" << s << "_I" << I << "_sig" << sigv << "_Ug" << U_g << "_Uc" << U_c << "_Udro" << Ut_dro << "_Ucel" << Ut_cel << "_Umam" << Ut_mam << "_Rg" << Rg << "_Rc" << Rc << "_Rep" << Rep << ".txt"; // results file naming convention
    nomF >> nomFichier; // Writing "nomF" to nomFichier file
    ofstream fout;
    fout.open(nomFichier);

	// population of N haploid individuals:
	ind * pop = new ind [Nv]; // "pop" chr pointer refers to a chr type (see "mutation.h") with Nv dimensions; "chr" has "gene", "cis", "trans", "sex" array
	ind * temp = new ind [Nv]; // temp for saving current generation during recombination
    ind * pc; // for recombination

    // array for measures from population: for each locus: sbarX, sbarY, hY, hXinact, attractXmale, attractXAct, attractXInact, attractY, dosY, dosXmale, dosXAct, dosXInact

    double ** measures;
    double * popAverages;
    if (output == 0)
    {
        measures = new double *[nbSv];
        for(i = 0; i < nbSv; i++)
            measures[i] = new double[8];
        popAverages = new double [5];
    }
    else
    {
        popAverages = new double [17];
    }

	// fitnesses:

	double * Wtot = new double [Nv]; // creates a new double pointer "Wtot" with dimension Nv -- each individual has a fitness!

	// for time length measure:

	time_t debut, fin; // "time_t" declares a time-calendar type; Debut-start, Fin-End
	struct tm *ptr; // sets "ptr" as a * for tm
	debut = time(0); // sets debut time to calendar time (not important)

    // initialization:

	for (i = 0; i < Nv; i++) // Loop through [Nv] times iterating over i (over each individual in the population)
    {
        // GENES
        pop[i].gene = new double [2*nbSv]; // two chromosomes
        temp[i].gene = new double [2*nbSv];

        // CIS regulators
        pop[i].cis = new double [2*nbSv]; // two chromosomes
        temp[i].cis = new double [2*nbSv];

        // TRANS regulators
        pop[i].trans = new double [2*3]; // two chromosomes, three trans reg's (m,n,o) corresponding to different strategies
        temp[i].trans = new double [2*3];

        for (k = 0; k < nbSv; k++)
        {
            pop[i].gene[k] = 1.0; // initialize all genes with 1 (WT/good) allele for Chromosome #1
            pop[i].gene[nbSv+k] = 1.0; // same for Chromosome #2
            pop[i].cis[k] = 1; // initialize all cis regulators with multiplier 0 for Chromosome #1 since exp(0)=1
            pop[i].cis[nbSv+k] = 1; // same for Chromosome #2
        }

        // Initialize trans multiplier coefficients: m = 0, n = 0, o = 0 since exp(0) = 1
        pop[i].trans[0] = 1; // Chrom 1
        pop[i].trans[1] = 1;
        pop[i].trans[2] = 1;

        pop[i].trans[3] = 1; // Chrom 2
        pop[i].trans[4] = 1;
        pop[i].trans[5] = 1;
    }

  	// generations:
	for (gen = 0; gen <= NbGenv; gen++) // Iterates over all generations
	{
        // Draw mutations:
        
        // Gene mutations
        mut = int(poisdev(UgTot)); // 2NuL number of expected mutations in total population

        for (nm = 0; nm < mut; nm++) // iterates over number of muts
        {
            indiv = rnd.randInt(N_1); // selects random individual; "rnd" in "MersenneTwister.h"
            site = rnd.randInt(twonbSm1); // selects random loci

            // draw mutational effect on fitness from exp distribution
            dg = s * (-log(rnd.rand())); // DISTRIBUTIONAL: Draw from exp dist of mean s, lambda = 1/s
//            dg = s; // CONSTANT: -s fitness mutation

            // "if" condition for biallele
           // if (pop[ind].gene[chrom*nbSv+site] == 1) // prevents cumulative mutations (i.e. biallelic constraint)
            //{
            //    pop[ind].gene[chrom*nbSv+site] += -dg; // always subtracts dg (away from optimal)
            //}

            // For a set s_max limit
            if (pop[indiv].gene[site] * (1 - dg) >= 1-s_max)
                pop[indiv].gene[site] *= (1 - dg);
            else
                pop[indiv].gene[site] = 1-s_max; // sets floor to 1-s_max
        }

        // Cis mutations
        mut = int(poisdev(UcTot)); // 2NuL number of expected mutations in total population
        for (nm = 0; nm < mut; nm++)
        {
            indiv = rnd.randInt(N_1);
            site = rnd.randInt(twonbSm1);

            // draw mutational effect on fitness from Gaussian distribution
            dc = sigv * gasdev(); // TURN THIS ON FOR: Continuous Mutation Effect
//            dc = s; // TURN THIS ON FOR: Fixed Mutation Effect
            pop[indiv].cis[site] += dc; // TURN THIS BACK ON!
        }

        // Trans mutations
        if (UTot_dro > 0)
        {
            mut = int(poisdev(UTot_dro)); // 2Nu number of expected mutations in total population
//           cout << "\n trans mut: " << mut;
            // ??? HOW TO DO ???
            for (nm = 0; nm < mut; nm++)
            {
                indiv = rnd.randInt(N_1);
                chrom = rnd.randInt(1);
 //              cout << "\n ind: " << ind << "  site: " << site << "  chrom: " << chrom;

                // draw mutational effect on fitness from Gaussian distribution
                dt = sigv * gasdev(); // TURN THIS ON FOR: Continuous Mutation Effect
                // dt = s; // TURN THIS ON FOR: Fixed Mutation Effect
                pop[indiv].trans[chrom*3+0] += dt; // mutate "m"
            }
        }
        if (UTot_cel > 0)
        {
            mut = int(poisdev(UTot_cel)); // 2Nu number of expected mutations in total population
            //           cout << "\n trans mut: " << mut;
            // ??? HOW TO DO ???
            for (nm = 0; nm < mut; nm++)
            {
                indiv = rnd.randInt(N_1);
                chrom = rnd.randInt(1);
                //              cout << "\n ind: " << ind << "  site: " << site << "  chrom: " << chrom;

                // draw mutational effect on fitness from Gaussian distribution
                dt = sigv * gasdev(); // TURN THIS ON FOR: Continuous Mutation Effect
                // dt = s; // TURN THIS ON FOR: Fixed Mutation Effect
                pop[indiv].trans[chrom*3+1] += dt; // mutate "m"
            }
        }
        if (UTot_mam > 0)
        {
            mut = int(poisdev(UTot_mam)); // 2Nu number of expected mutations in total population
            //           cout << "\n trans mut: " << mut;
            // ??? HOW TO DO ???
            for (nm = 0; nm < mut; nm++)
            {
                indiv = rnd.randInt(N_1);
                chrom = rnd.randInt(1);
                //              cout << "\n ind: " << ind << "  site: " << site << "  chrom: " << chrom;

                // draw mutational effect on fitness from Gaussian distribution
                dt = sigv * gasdev(); // TURN THIS ON FOR: Continuous Mutation Effect
                // dt = s; // TURN THIS ON FOR: Fixed Mutation Effect
                pop[indiv].trans[chrom*3+2] += dt; // mutate "m"
            }
        }

        wmmax = 0;
        for (i = 0; i < Nmales; i++) // Loops through each male in the population
        {
            w = Wmale(pop[i], expdom, s_max, Q0, I, nbSv);
            Wtot[i] = w;
            if (wmmax < w)
                wmmax = w;
        }

        wfmax = 0;
        for (i = Nmales; i < Nv; i++) // Loops through each female in the population
        {
            w = Wfemale(pop[i], expdom, s_max, Q0, I, nbSv);
            Wtot[i] = w;
            if (wfmax < w)
                wfmax = w;
        }

        Nmales_1 = Nmales-1;
        Nfemales_1 = Nv-Nmales-1;
        NjuvM = 0;
        NjuvF = 0;

        // measures phenotypic moments and writes in result file every "pasv" generations:


        if (gen % pasv == 0) // "gen" is the generational index (iterable) and once it reaches a multiple of the time interval per generation "pasv", write results
        {
            if (output == 0)
            {
                record_output(pop, Wtot, measures, popAverages, nbSv, Nmales, Nv, expdom);

                fout << gen << " " << popAverages[0] << " " << popAverages[1] << " " << popAverages[2] << " " << popAverages[3] << " " << popAverages[4] << " " ;
                for (j = 0; j < nbSv; j++)
                    for (i = 0; i < 8; i++)
                        fout << measures[j][i] << " ";
                fout << endl;
            }
            else
            {
                record_averages(pop, Wtot, popAverages, nbSv, Nmales, Nv, expdom, s_max);
                
                fout << gen << " ";
                for (i = 0; i < 17; i++)
                {
                    fout << popAverages[i] << " ";
                    allAverages[gen/pasv][i] = ((allAverages[gen/pasv][i])*Rep + popAverages[i]) / (Rep+1);
                }
                fout << endl;
            }
        }


        // next generation (meiosis):
		for (j = 0; j < Nv; j++)
		{
            do{
				p1 = rnd.randInt(Nmales_1); // select a father

			} while (rnd.rand() > (Wtot[p1] / wmmax));
            // second parent:

			do{
                p2 = rnd.randInt(Nfemales_1); // select random parent 2

            } while (rnd.rand() > (Wtot[Nmales + p2] / wfmax));

            // Sex determination
            if (rnd.rand() < 0.5) // offspring is male
                off_sex = 0;
            else
                off_sex = 1;

            // recombination:

            if (gen < NbPrelimv)
            {
                if (off_sex == 0)
                {
                    recInit(temp[NjuvM], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvM += 1;
                }
                else
                {
                    recInit(temp[Nv-1-NjuvF], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvF += 1;
                }
            }
            else
            {
                if (off_sex == 0)
                {
                    rec(temp[NjuvM], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvM += 1;
                }
                else
                {
                    rec(temp[Nv-1-NjuvF], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvF += 1;
                }
            }
		}

        // update population:

        pc = pop;
        pop = temp;
        temp = pc;

        Nmales = NjuvM;
	}

    fin = time(0); // NOT IMPORTANT (C function for getting the calendar time to write to the results file

    // writes in output file:

    fprintf(fichierS, "\n\nResults in file ");
    fprintf(fichierS, nomFichier);
    fprintf(fichierS, "\n");

    // time length:

    int temps = int(difftime(fin, debut));
    fprintf(fichierS, "\n%d generations took %d hour(s) %d minute(s) %d seconds\n", NbGenv, temps / 3600, (temps % 3600) / 60, temps % 60);

    // date and time:

    ptr=localtime(&fin); //  & = address of; localtime(&fin) gives the corresponding time zone time for "fin"
    fprintf(fichierS, asctime(ptr)); // ERIC: prints asctime(ptr) to the FILE stream with pointer fichierS; asctime(time) gives "time" in human-readable format

	// FREES MEMORY AFTER RUNNING
    for (i = 0; i < Nv; i++)
    {
        // delete [] "pointer" deletes/frees the array memory for the given pointer
        delete [] pop[i].gene;
        delete [] temp[i].gene;
        delete [] pop[i].trans;
        delete [] temp[i].trans;
        delete [] pop[i].cis;
        delete [] temp[i].cis;
    }
	delete [] pop;
    delete [] temp;
    delete [] Wtot;
    if (output == 0)
    {
        for(i = 0; i < nbSv; i++)
            delete [] measures[i];
        delete [] measures;
    }
    delete [] popAverages;
}
