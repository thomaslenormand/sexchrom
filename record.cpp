#include "mutation.h"
#include <cmath>
using namespace std;

void record_output(ind * pop, double * Wtot, double ** measures, double * popAve, int nbSv, int Nmales, int Nv, double expdom)
{
    // measures from population: for each locus: sbarX, sbarY, hY, hXinact, attractXmale, attractY, attractXAct, attractXInact

    double WbarMales, WbarFemales, mbar, nbar, obar, y, transDro, transCel, transMam, c1, c2;
    int i, j;
    int twoNfemale = 2*(Nv - Nmales);

    for (j = 0; j < nbSv; j++)
        for (i = 0; i < 8; i++)
            measures[j][i] = 0;

    WbarMales = 0;
    mbar = 0; nbar = 0; obar = 0;
    for (i = 0; i < Nmales; i++)
    {
        WbarMales += Wtot[i];
        
        transDro = (((pop[i].trans[0] > 0) ? pop[i].trans[0] : 0) + ((pop[i].trans[3] > 0) ? pop[i].trans[3] : 0))/2.0;
        
        mbar += transDro;
        nbar += (((pop[i].trans[1] > 0) ? pop[i].trans[1] : 0) + ((pop[i].trans[4] > 0) ? pop[i].trans[4] : 0))/2.0;
        obar += (((pop[i].trans[2] > 0) ? pop[i].trans[2] : 0) + ((pop[i].trans[5] > 0) ? pop[i].trans[5] : 0))/2.0;
        
        for (j = 0; j < nbSv; j++)
        {
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);
            
            measures[j][0] += 1 - pop[i].gene[nbSv+j];
            measures[j][1] += 1 - pop[i].gene[j];
            y = c1 / (c1 + c2*transDro);
            measures[j][2] += pow(y, expdom);
            measures[j][4] += c2*transDro;
            measures[j][5] += c1;
        }
    }
    WbarMales /= Nmales;

    WbarFemales = 0;
    for (i = Nmales; i < Nv; i++)
    {
        WbarFemales += Wtot[i];
        
        transMam = (((pop[i].trans[2] > 0) ? pop[i].trans[2] : 0) + ((pop[i].trans[5] > 0) ? pop[i].trans[5] : 0))/2.0;
        transCel = (((pop[i].trans[1] > 0) ? pop[i].trans[1] : 0) + ((pop[i].trans[4] > 0) ? pop[i].trans[4] : 0))/2.0;
        
        mbar += (((pop[i].trans[0] > 0) ? pop[i].trans[0] : 0) + ((pop[i].trans[3] > 0) ? pop[i].trans[3] : 0))/2.0;
        nbar += transCel;
        obar += transMam;

        for (j = 0; j < nbSv; j++)
        {
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);
            
            measures[j][0] += 2 - pop[i].gene[j] - pop[i].gene[nbSv+j];
            y = (c1*transMam) / (c1*transMam + c2);
            measures[j][3] += pow(y, expdom);
            y = (c2*transMam) / (c2*transMam + c1);
            measures[j][3] += pow(y, expdom);
            measures[j][6] += (c1 + c2) * transCel;
            measures[j][7] += (c1 + c2) * transMam * transCel;
        }

    }
    for (j = 0; j < nbSv; j++)
    {
        measures[j][0] /= (Nmales + twoNfemale); // sbarX
        measures[j][1] /= Nmales; // sbarY
        measures[j][2] /= Nmales; // hY
        measures[j][3] /= twoNfemale; // hXInact
        measures[j][4] /= Nmales; // attractXMales
        measures[j][5] /= Nmales; // attractY
        measures[j][6] /= twoNfemale; // attractXAct
        measures[j][7] /= twoNfemale; // attractXInact
    }
    WbarFemales /= (Nv - Nmales);
    mbar /= Nv;
    nbar /= Nv;
    obar /= Nv;

    popAve[0] = WbarMales;
    popAve[1] = WbarFemales;
    popAve[2] = mbar;
    popAve[3] = nbar;
    popAve[4] = obar;
}

void record_averages(ind * pop, double * Wtot, double * popAve, int nbSv, int Nmales, int Nv, double expdom, double smax)
{
    // measures from population: Wbarmale, WbarFemale, mbar, nbar, obar, sbarX, sbarY, hY, hXinact, attractXmale, attractY, attractXAct, attractXInact,
    // number of loci half-dead, dead, half-silenced, silenced on Y
    
    double WbarMales, WbarFemales, mbar, nbar, obar, y, transDro, transCel, transMam, c1, c2;
    int i, j;
    int twoNfemale = 2*(Nv - Nmales);
    double halfminW = (1-smax/2.0);
    double minW = 1-smax;
    int cmptyM = 0;
    int cmptyF = 0;
    
    for (i = 0; i < 17; i++)
        popAve[i] = 0;
    
    WbarMales = 0;
    mbar = 0; nbar = 0; obar = 0;
    for (i = 0; i < Nmales; i++)
    {
        WbarMales += Wtot[i];
        
        transDro = (((pop[i].trans[0] > 0) ? pop[i].trans[0] : 0) + ((pop[i].trans[3] > 0) ? pop[i].trans[3] : 0))/2.0;
        
        mbar += transDro;
        nbar += (((pop[i].trans[1] > 0) ? pop[i].trans[1] : 0) + ((pop[i].trans[4] > 0) ? pop[i].trans[4] : 0))/2.0;
        obar += (((pop[i].trans[2] > 0) ? pop[i].trans[2] : 0) + ((pop[i].trans[5] > 0) ? pop[i].trans[5] : 0))/2.0;
        
        for (j = 0; j < nbSv; j++)
        {
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);
            
            popAve[5] += 1 - pop[i].gene[nbSv+j];
            popAve[6] += 1 - pop[i].gene[j];
            if (c1 + c2*transDro > 0)
            {
                y = c1 / (c1 + c2*transDro);
                popAve[7] += pow(y, expdom);
                cmptyM++;
            }
            popAve[9] += c2*transDro;
            popAve[10] += c1;
            if (pop[i].gene[j] < halfminW)
                popAve[13] += 1;
            if (pop[i].gene[j] == minW)
                popAve[14] += 1;
            if (y < 0.25)
                popAve[15] += 1;
            if (y < 0.01)
                popAve[16] += 1;
        }
    }
    WbarMales /= Nmales;
    
    WbarFemales = 0;
    for (i = Nmales; i < Nv; i++)
    {
        WbarFemales += Wtot[i];
        
        transMam = (((pop[i].trans[2] > 0) ? pop[i].trans[2] : 0) + ((pop[i].trans[5] > 0) ? pop[i].trans[5] : 0))/2.0;
        transCel = (((pop[i].trans[1] > 0) ? pop[i].trans[1] : 0) + ((pop[i].trans[4] > 0) ? pop[i].trans[4] : 0))/2.0;
        
        mbar += (((pop[i].trans[0] > 0) ? pop[i].trans[0] : 0) + ((pop[i].trans[3] > 0) ? pop[i].trans[3] : 0))/2.0;
        nbar += transCel;
        obar += transMam;
        
        for (j = 0; j < nbSv; j++)
        {
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);
            
            popAve[5] += 2 - pop[i].gene[j] - pop[i].gene[nbSv+j];
            if (c1*transMam + c2 > 0)
            {
                y = (c1*transMam) / (c1*transMam + c2);
                popAve[8] += pow(y, expdom);
                cmptyF++;
            }
            if (c2*transMam + c1 > 0)
            {
                y = (c2*transMam) / (c2*transMam + c1);
                popAve[8] += pow(y, expdom);
                cmptyF++;
            }
            popAve[11] += (c1 + c2) * transCel;
            popAve[12] += (c1 + c2) * transMam * transCel;
        }
        
    }
    popAve[5] /= (nbSv *(Nmales + twoNfemale)); // sbarX
    popAve[6] /= (nbSv * Nmales); // sbarY
    popAve[7] /= cmptyM; // hY
    popAve[8] /= cmptyF; // hXInact
    popAve[9] /= (nbSv * Nmales); // attractXMales
    popAve[10] /= (nbSv * Nmales); // attractY
    popAve[11] /= (nbSv * twoNfemale); // attractXAct
    popAve[12] /= (nbSv * twoNfemale); // attractXInact
    popAve[13] /= (nbSv * Nmales);
    popAve[14] /= (nbSv * Nmales);
    popAve[15] /= (nbSv * Nmales);
    popAve[16] /= (nbSv * Nmales);
    
    WbarFemales /= (Nv - Nmales);
    mbar /= Nv;
    nbar /= Nv;
    obar /= Nv;
    
    popAve[0] = WbarMales;
    popAve[1] = WbarFemales;
    popAve[2] = mbar;
    popAve[3] = nbar;
    popAve[4] = obar;
}



