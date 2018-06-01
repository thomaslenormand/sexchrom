#include "mutation.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

// rec function: generates recombinant genome segment "res"
// from parental genome segments "c1" and "c2"
// nbCo is the number of cross-overs in the genome segment,
// nS the number of selected loci, nU nb of loci affecting U


void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex)
{
    vector<double> posCoG;
    int i, j, locus, nbCo, ChrInit, strand;

    // paternally inherited gamete
    double Rsex = Rg; // NEW linkage between the first cis reg and the sex locus. Only matters in male where the sex locus is heterozygous. Set to Rg to avoid adding a parameter
    nbCo = poisdev(ML+Rsex); // NEW number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    if (off_sex == 0) // if offspring is a male, it inherits the sex-determining locus from the proto-Y
        ChrInit = 0;
    else
        ChrInit = 1;

    // CIS REG:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while (locus * Rg + Rsex < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        {
            offspring.cis[locus] = father.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    {
        offspring.cis[locus] = father.cis[strand + locus];
        locus++;
    }

    // GENES:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus * Rg + Rc + Rsex) < posCoG[j]) // NEW note that position of gene i is i*Rg + Rc
        {
            offspring.gene[locus] = father.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[locus] = father.gene[strand + locus];
        locus++;
    }

    // maternally inherited gamete

    nbCo = poisdev(ML); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    ChrInit = rnd.randInt(1); // which X chromosome contributes at position 0

    // CIS REG:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while (locus * Rg < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        {
            offspring.cis[nS+locus] = mother.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    {
        offspring.cis[nS+locus] = mother.cis[strand + locus];
        locus++;
    }

    // GENES:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus * Rg + Rc) < posCoG[j]) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[nS+locus] = mother.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[nS+locus] = mother.gene[strand + locus];
        locus++;
    }

    // TRANS REGULATORS

    for (i = 0; i < 3; i++)
    {
        // paternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i] = father.trans[ChrInit*3+i];

        // maternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i+3] = mother.trans[ChrInit*3+i];
    }
}



void rec(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex)
{
    vector<double> posCoG;
    int i, j, locus, nbCo, ChrInit, strand;

    if (off_sex == 0) // if offspring is a male
    {
        for (i = 0; i < nS; i++)
        {
            offspring.gene[i] = father.gene[i]; // father transmits his Y chromosome
            offspring.cis[i] = father.cis[i];
        }
    }
    else
    {
        for (i = 0; i < nS; i++)
        {
            offspring.gene[i] = father.gene[nS+i]; // father transmits his X chromosome
            offspring.cis[i] = father.cis[nS+i];
        }
    }

    nbCo = poisdev(ML); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end()); // ERIC: sort indices for crossovers -- first to last

    ChrInit = rnd.randInt(1); // which X chromosome contributes at position 0

    // CIS REG:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while (locus * Rg < posCoG[j]) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
        {
            offspring.cis[nS+locus] = mother.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
    {
        offspring.cis[nS+locus] = mother.cis[strand + locus];
        locus++;
    }

    // GENES:
    locus = 0; // ERIC: initialize locus again
    for (j = 0; j < nbCo; j++) // ERIC: iterate over the number of crossovers
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus * Rg + Rc) < posCoG[j]) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[nS+locus] = mother.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[nS+locus] = mother.gene[strand + locus];
        locus++;
    }

    // TRANS REGULATORS

    for (i = 0; i < 3; i++)
    {
        // paternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i] = father.trans[ChrInit*3+i];

        // maternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.trans[i+3] = mother.trans[ChrInit*3+i];
    }
}
