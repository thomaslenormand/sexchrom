// Functions to open input and output files,
// read parameter values from input file and
// write them in output file.

#include "mutation.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;

// opens input file:

void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r"); // ERIC: Opens fichierLecture file in read mode
}


// opens output file:

void ouvrirFichierS()
{
	fichierS = fopen(fichierEcriture,"a"); // ERIC: Opens fichierEcriture file in append mode
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0
bool lireFichier(int &Nr, double &sigr, int &nbSr, int &NbGenr, int &NbPrelimr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &Ut_dro, double &Ut_cel, double &Ut_mam, double &Rgr, double &Rcr, int &Rep, int &outputr)
{
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
	{
		term = true;
	}
	else
	{
		fscanf(fichierE,"%d ",&Nr);  // ERIC: "fscanf(FILE_ptr, format, address)" reads data from stream and stores at address in given format
		fscanf(fichierE,"%lf ",&sigr);
		fscanf(fichierE,"%d ",&nbSr);
		fscanf(fichierE,"%d ",&NbGenr);
        fscanf(fichierE,"%d ",&NbPrelimr);
		fscanf(fichierE,"%d ",&pasr);
        fscanf(fichierE,"%lf ",&sr);
        fscanf(fichierE,"%lf ",&s_maxr);
        fscanf(fichierE,"%lf ",&Ir);
        fscanf(fichierE,"%lf ",&U_gr);
        fscanf(fichierE,"%lf ",&U_cr);
        fscanf(fichierE,"%lf ",&Ut_dro);
        fscanf(fichierE,"%lf ",&Ut_cel);
        fscanf(fichierE,"%lf ",&Ut_mam);
        fscanf(fichierE,"%lf ",&Rgr);
        fscanf(fichierE,"%lf ",&Rcr);
        fscanf(fichierE,"%d ",&Rep);
        fscanf(fichierE,"%d ",&outputr);

		term = false;
	}
	return term;
}


// writes parameter values in output file:

void ecrireParametres(int Nv, double sigv, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double Ut_dro, double Ut_cel, double Ut_mam, double Rgv, double Rcv, int outputv)
{
	fprintf(fichierS,"\n_________________________________________\n"); // ERIC: "fprintf(FILE_ptr, string_format, var)" writes (string_format, var) to file stream
	fprintf(fichierS,"\nN = %d", Nv);
	fprintf(fichierS,", sigma = %g", sigv);
	fprintf(fichierS,", nbS = %d", nbSv);
	fprintf(fichierS,"\ngenerations = %d", NbGenv);
    fprintf(fichierS,", NbPrelim = %d", NbPrelimv);
	fprintf(fichierS,", pas = %d", pasv);
    fprintf(fichierS,", s = %g", s);
    fprintf(fichierS,", s_max = %g", s_max);
    fprintf(fichierS,", s = %g", I);
    fprintf(fichierS,", U_g = %g", U_g);
    fprintf(fichierS,", U_c = %g", U_c);
    fprintf(fichierS,", Ut_dro = %g", Ut_dro);
    fprintf(fichierS,", Ut_cel = %g", Ut_cel);
    fprintf(fichierS,", Ut_mam = %g", Ut_mam);
    fprintf(fichierS,", Rg = %g", Rgv);
    fprintf(fichierS,", Rc = %g", Rcv);
    fprintf(fichierS,", output = %d", outputv);
}
