/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: PNRdyn.cpp
/
/ Function: PNRdyn()
/
/ Summary:
/
/ This function computes the probability of ruin between the two RF(t) buckets provided using a dynamic program. The function accepts an array of
/ probabilities from the prior time-point along with a vector of bucket numbers that point to the buckets that represent unique probabilities (the last
/ bucket in the sequence). The probabilities computed are entered into the array passed to this function as a pointer. If this array already holds
/ values they will be replaced by the new computations, it is not required that the array be empty. The array of prior probabilities is named Vp[] and
/ the array of current time-point probabilities is V[]. To compute probabilities using a dynamic program only the CDF of the random variables
/ representing inflation/expense-adjusted returns is needed. When computing a standard ruin probability the inflation/expense-adjusted returns are
/ assumed normally distributed for this implementation with means/variances/covariances of the corresponding stock/bond returns specified in the control
/ file. When computing gradient entries the special density g() with CDF G() is used to represent the inflation/expense-adjusted return at the current
/ time-point. When computing diagonal Hessian entries the densities h1() and h2() with CDFs H1() and H2() are used to represent inflation/expense-
/ adjusted returns. Off-diagonal Hessian entries use the gradient special density g() with CDF G() at both time-points represented by the diagonal
/ indices.
/
/ Parameters:
/
/ 1.) Array of precomputed moments & constants for the given time-point and equity ratio (to save processing time): m, mp, v, vp, vpp, mva, kh1.
/ 2.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/ 3.) A 3-value array holding the start ruin factor bucket, end ruin factor bucket, and total # buckets in the discretization. The probability of ruin
/ is computed and entered into the incoming array for all buckets between the start and end buckets (inclusive). Allowing for bucket-specific calls
/ enables multi-threading when this function is invoked from a wrapper.
/ 4.) An array (as pointer) of long doubles holding the probabilities of ruin for each bucket at the prior time-point. The size of this array is the
/ total number of buckets in the discretization (i.e., RFMax*Prec).
/ 5.) An array (as pointer) of long doubles that will be populated between the buckets given in parameter #3. The size of this array is the total number
/ of buckets in the discretization (i.e., RFMax*Prec).
/ 6.) A vector containing the bucket numbers from the prior time-point where unique ruin probabilities reside. If there is a string of consecutive
/ buckets with the same ruin probabilities they can be treated as a single larger bucket to save processing time at the current time-point.
/ 7.) Vector of constants. All CDFs in this application can be expressed as a linear combination of known CDF calls. This is the constant vector for
/ that linear combination. This vector along with the vector of gamma random variables in the next parameter are passed to the GetCDF() function.
/ 8.) A vector of gamma distributed random variables (class type=boost distribution) needed to build the CDF for a given time-point and equity ratio.
/ The CDF G() of g(), and the CDF H1() of h1() both use 2 gamma RVs. The CDF H2() of h2() uses 4 gamma RVs. This vector along with the constant
/ vector from the previous argument are passed to the GetCDF() function.
/
/ Return Value:
/
/ This function returns no value but populates the array supplied in the 5th argument with the corresponding ruin probabilities between the two
/ RF(t) buckets specified in the 3rd argument.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Collections.Generic;
using Accord.Statistics.Distributions.Univariate;

namespace OptimalEquity.Model
{
    public static class PNRdyn
    {
        public static void Run(double[] mts, int prec, int[] bkts, double[] vp, double[] v, List<int> prB,
            List<double> cLocal, List<GammaDistribution> g)
        {
            const double tiethresh = 0.5;
            int nuqbkts = prB.Count;
            int ties = 0;
            int cont = 1;

            NormalDistribution normdist = new NormalDistribution(mts[0], Math.Sqrt(mts[2]));

            // Iterate over buckets until a probability of 1.00 is encountered.
            for (int b = bkts[0]; cont == 1 && b <= bkts[1]; ++b)
            {
                double rf = (double) b / prec;
                double pruin;

                // Derive P(Ruin) for this bucket, time point, and asset allocation.
                double cdfval = GetCDF.Run(mts, cLocal, normdist, g, rf);
                double eprob;
                if (cdfval >= 1.00)
                {
                    eprob = vp[bkts[2] - 1];
                }
                else
                {
                    var rhsCdf = 1.00;
                    double lhsCdf = GetCDF.Run(mts, cLocal, normdist, g, rf * (1 + (prec / 1.5)));
                    eprob = (rhsCdf - lhsCdf) * vp[0]; // First bucket, unique processing.
                    rhsCdf = lhsCdf;
                    for (int pb = 2;
                        pb <= nuqbkts;
                        ++pb) // All others but last bucket, standard processing for unique probs only.
                    {
                        lhsCdf = GetCDF.Run(mts, cLocal, normdist, g, rf * (1.0 + prec / (prB[pb - 1] + 0.5)));
                        eprob = eprob + (rhsCdf - lhsCdf) * vp[prB[pb - 1] - 1];
                        rhsCdf = lhsCdf;
                    }
                    eprob = (eprob + (rhsCdf - cdfval) * vp[bkts[2] - 1]) /
                            (1.00 - cdfval); // Last bucket, unique processing, make it conditional.
                }

                // Deal with numerical instability near zero and one.
                if (ties == 0)
                {
                    pruin = cdfval + eprob - cdfval * eprob;
                    if (pruin > tiethresh)
                    {
                        pruin = 1.00 - (1.00 - cdfval) * (1.00 - eprob);
                        ties = 1;
                    }
                }
                else
                {
                    pruin = 1.00 - (1.00 - cdfval) * (1.00 - eprob);
                }

                // Load ruin probability into array.
                v[b - 1] = pruin;

                // Pruning: Once we hit probability of 1.00, set the remaining buckets manually.
                if (v[b - 1] >= 1.00)
                {
                    for (int bb = b + 1; bb <= bkts[1]; ++bb)
                        v[bb - 1] = 1.00;
                    cont = 0;
                }
            }
        }
    }
}
