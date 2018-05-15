#include "mutation.h"
#include <cmath>
using namespace std;

extern MTRand rnd;

double Wmale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j, allele1, allele2;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus;

    double transDro = (((parent.trans[0] > 0) ? parent.trans[0] : 0) + ((parent.trans[3] > 0) ? parent.trans[3] : 0))/2; // the effect of each allele is set to zero if negative

    for (j = 0; j < nbS; j++)
    {
        allele1 = j;
        allele2 = nbS+j;

        e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);  // cis_y
        e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0) * transDro;  // (m1+m2)/2 * cis_x

        explv = e1 + e2; // (m1+m2)/2 * cis_x + cis_y
        //stabilizing selection fitness effect
        //                   w = w * exp(-I*pow((explv-Qopt),2)); // OLD

        if (explv == 0)
            w *= 1 - smax;

        else
        {
            //effectLocus = 1-smax*(1-exp(-I*pow((explv-Qopt),2))); // OLD
            effectLocus = 1-smax*(1-explv*exp(1-explv/Qopt)/Qopt);      // NEW

            if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e1/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[nbS+j] + h * (parent.gene[j] - parent.gene[nbS+j]); <-- additive fitness
                effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2]));
            }
            else // if first allele is fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e2/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[j] + h * (parent.gene[nbS+j] - parent.gene[j]); <-- additive fitness
                effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1]));
            }
            if (effectLocus < 1 - smax)
                w *= 1 - smax;
            else
                w *= effectLocus;
        }
    }
    return w;
}


double Wfemale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus;
    int allele1, allele2;

    double transMam = (((parent.trans[2] > 0) ? parent.trans[2] : 0) + ((parent.trans[5] > 0) ? parent.trans[5] : 0))/2;
    double transCel = (((parent.trans[1] > 0) ? parent.trans[1] : 0) + ((parent.trans[4] > 0) ? parent.trans[4] : 0))/2;

    int chrom = rnd.randInt(1);// pick random x

    if (chrom == 0) // if first chromosome randomly selected
    {
        for (j = 0; j < nbS; j++)
        {
            allele1 = j;
            allele2 = nbS+j; // to avoid doing the sum multiple times

            // (n1+n2)/2 * [(o1+o2)/2*cis_x1 + cis_x2]
            e1 = transMam * ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);
            e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0);
            explv = (e1 + e2) * transCel;

            if (explv == 0)
                w *= 1 - smax;

            else
            {
                //stabilizing selection fitness effect
                //                    w = w * exp(-I*pow((explv-Qopt),2)); // OLD
                //effectLocus = 1-smax*(1-exp(-I*pow((explv-Qopt),2))); // OLD
                effectLocus = 1-smax*(1-explv*exp(1-explv/2)/2);      // NEW

                if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e1/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2])); // product fitnesses
                }
                else // if first allele fittest
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e2/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1])); // product fitnesses
                }
                if (effectLocus < 1 - smax)
                    w *= 1 - smax;
                else
                    w *= effectLocus;
            }
        }
    }

    else
    {
        for (j = 0; j < nbS; j++)
        {
            allele1 = j;
            allele2 = nbS+j; // to avoid doing the sum multiple times

            // (n1+n2)/2 * [(o1+o2)/2*cis_x2 + cis_x1]
            e2 = transMam * ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0);
            e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);
            explv = (e1 + e2) * transCel;

            if (explv == 0)
                w *= 1 - smax;

            else
            {
                //stabilizing selection fitness effect
                //                    w = w * exp(-I*pow((explv-Qopt),2)); // OLD
                //effectLocus = 1-smax*(1-exp(-I*pow((explv-Qopt),2))); // OLD
                effectLocus = 1-smax*(1-explv*exp(1-explv/2)/2);      // NEW

                if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
                    // (n1+n2)/2 * [(o1+o2)/2*cis_x2 + cis_x1]
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e1/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2])); // product fitnesses
                }
                else // if first allele fittest
                {
                    // h is the dominance coefficient, H0 = 0.25
                    //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                    h = pow((e2/(e1+e2)), expdom);

                    effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1])); // product fitnesses
                }
                if (effectLocus < 1 - smax)
                    w *= 1 - smax;
                else
                    w *= effectLocus;
            }
        }
    }
    return w;
}



