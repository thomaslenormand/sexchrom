// Header file: definitions of global variables, function prototypes

#ifndef MUTATION_H
#define MUTATION_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// Global variables:

#define fichierLecture "parametres.txt"     // names of input
#define fichierEcriture "resultats.txt"		// and output files

// "chr": represents a chromosome

struct ind
{
    double * gene; // allelic values at loci corresponding to fitness
    double * cis; // cistronic values at loci controlling regulation [doubles]
    double * trans; // trans values [m, n, o] for drosophila, c elegans, mammal ()
    //   	Q_male = (m1+m2)/2 * q_x + q_y
    //      Q_female = [(o1+o2)/2 * q1 + q2] * (n1+n2)/2

};

// Function prototypes:

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(int Nv, double sigv, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double Ut_all, double Ut_male, double Ut_dro, double Ut_cel, double Ut_mam, double Rgv, double Rgc, int outputv);
bool lireFichier(int &Nr, double &sigr, int &nbSr, int &NbGenr, int &NbPrelimr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &Ut_all, double &Ut_male, double &Ut_dro, double &Ut_cel, double &Ut_mam, double &Rgr, double &Rcr, int &Rep, int &outputr);
void recursion(int Nv, double sigv, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double Ut_all, double Ut_male, double Ut_dro, double Ut_cel, double Ut_mam, double Rg, double Rc, int Rep, int output, double** allAverages);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex);
void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex);
double Wmale(ind &parent, double h0, double s_max, double Q0, double I, int nbS);
double Wfemale(ind &parent, double h0, double s_max, double Q0, double I, int nbS);
void record_output(ind * pop, double * Wtot, double ** measures, double * popAve, int nbSv, int Nmales, int Nv, double h0);
void record_averages(ind * pop, double * Wtot, double * popAve, int nbSv, int Nmales, int Nv, double expdom, double smax);

#endif
