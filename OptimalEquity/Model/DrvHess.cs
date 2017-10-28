/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: DrvHess.cpp
/
/ Function: DrvHess()
/
/ Summary:
/
/ This function computes the Hessian values and populates, then returns, a square TDxTD symmetric matrix with these values. Each element of the Hessian
/ matrix can be expressed as a linear combination of ruin probabilities. The form of this linear combination differs if the element is on the diagonal
/ or not. The user selects an estimation method of simulation (sim) or dynamic programming (dp) for estimating these ruin probabilities. The off-
/ diagonal elements use first order partial derivatives in both the i,j dimension and the diagonal elements use 2nd order partial derivatives in the i-i
/ dimension. The off-diagonal elements can make use of the existing gradient vector and success probability therefore only requires one new ruin
/ probability be estimated which uses the gradient special density g() in both the i and j dimensions. The diagonal elements can make use of the success
/ probability only and requires 2 new ruin probabilities be estimated with either simulation or a dynamic program. The diagonal elements make use of the
/ special densities h1() and h2(), or their corresponding CDFs.
/
/ Parameters:
/
/ 1.) Estimation type: Either dynamic programming (dp) or simulation (sim).
/ 2.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 3.) Array of long doubles holding the glide-path to use for this calculation (as pointer). The corresponding array must have exactly TD elements.
/ 4.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 5.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 6.) The sample size to use when simulating the ruin/success probability. (Note that the ruin probability is the number of times ruin occurs divided
/ by the number of trial runs. The success probability is 1 minus the ruin probability. (Here the probability will use the special densities h1()
/ and h2() when building diagonal Hessian elements and will use the gradient density g() for off-diagonal elements.) This parameter is only used
/ when estimation type is simulation.
/ 7.) The number of buckets to use when discretizing the ruin factor dimension. This value equals RFMax*(Discretization Precision), where these 2 values
/ are set by the user in the control file. (This parameter is only used when type="dp".)
/ 8.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/ 9.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/ 10.) A pointer to the gradient vector for this glide-path. The off-diagonal Hessian elements include ruin calculations that can be derived directly
/ from the existing gradient values and these are used instead of recalculating them.
/ 11.) The success probability for the given glide-path. This value is also used for computing off-diagonal Hessian elements and is not recomputed each
/ time.
/
/ Return Value:
/
/ This function returns a square TDxTD symmetric matrix representing the Hessian of the glide-path supplied.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class DrvHess
    {
        public static double[,] Run(string type, double[] prms, double[] a, int fxTD, double rf0, long n, int nbuckets, int prec, int plproc, double[] gradvctr, double stdpnr)
        {
            double[,] hess = new double[fxTD, fxTD];

            // Derive the Hessian matrix.
            Trace.WriteLine("");
            Trace.Write("Building Hessian ");
            for (int i = 0; i < fxTD; ++i)
            {
                var blnks = new String(' ', (i > 0 ? 17 : 0) + i);
                Trace.Write(blnks);
                for (int j = 0; j < fxTD; ++j)
                {
                    if (j >= i)
                    {
                        Trace.Write(".");
                        if (i < 0 || j < 0 || i >= fxTD || j >= fxTD)
                        {
                            Trace.Write("ERROR: Both i and j must be integers between 0 and ");
                            Trace.Write(fxTD - 1);
                            Trace.Write(", i=");
                            Trace.Write(i);
                            Trace.Write(" and j=");
                            Trace.WriteLine(j);
                            Trace.WriteLine("EXITING...DrvHess()...");
                            Console.Read();
                            Environment.Exit(1);
                        }
                        else if (i != j) // ***** Off-diagonal elements ***** //
                        {
                            // Define needed quantities and return the off-diagonal Hessian element.
                            //========================================================================
                            double k1 = Funcs.vp(prms, a[i]) / (2.00 * Funcs.v(prms, a[i])) +
                                        Math.Pow(Funcs.mp(prms), 2) / (2.00 * Funcs.vp(prms, a[i]));
                            double k2 = Funcs.vp(prms, a[j]) / (2.00 * Funcs.v(prms, a[j])) +
                                        Math.Pow(Funcs.mp(prms), 2) / (2.00 * Funcs.vp(prms, a[j]));
                            int[] prtls = {i, j, -1, -1};
                            hess[i, j] =
                                k1 * k2 * (GetPNR.Run(type, prms, a, fxTD, rf0, n, nbuckets, prec, prtls, plproc) -
                                           gradvctr[i] / k1 - gradvctr[j] / k2 - stdpnr);
                        }
                        else // ***** Diagonal elements ***** //
                        {
                            // Define needed quantities and return the diagonal Hessian element.
                            // [Note #1: K1 is for h1, K2 is for h2, and K3 is for f(). And the diagonal term is K1*h1() + K2*h2() + K3*f().]
                            // [Note #2: The Hessian diagonals will divide by zero for one alpha, check for that value and exit if encountered.]
                            if (Math.Abs(Funcs.v(prms, a[i]) * Funcs.vpp(prms) - 2.00 * Math.Pow(Funcs.vp(prms, a[i]), 2)) < 1e-15)
                            {
                                Trace.Write("ERROR: Hessian diagonal element does not exist for alpha value=");
                                Trace.WriteLine($"{a[i]} encountered at time point t={i}.");
                                Trace.WriteLine("");
                                Trace.WriteLine("EXITING...DrvHess()...");
                                Console.Read();
                                Environment.Exit(1);
                            }
                            else
                            {
                                double k1 = (Funcs.v(prms, a[i]) + Math.Pow(Funcs.kh1(prms, a[i]), 2)) *
                                            (Funcs.v(prms, a[i]) * Funcs.vpp(prms) - 2.00 * Math.Pow(Funcs.vp(prms, a[i]), 2)) /
                                            (2.00 * Math.Pow(Funcs.v(prms, a[i]), 3));
                                double k2 = (Math.Pow(Funcs.vp(prms, a[i]), 2) + 2.00 * Funcs.v(prms, a[i]) * Math.Pow(Funcs.mp(prms), 2)) /
                                            (2.00 * Math.Pow(Funcs.v(prms, a[i]), 2));
                                double k3 = -((Funcs.vpp(prms) * Funcs.v(prms, a[i]) - Math.Pow(Funcs.vp(prms, a[i]), 2) +
                                               2.00 * Funcs.v(prms, a[i]) * Math.Pow(Funcs.mp(prms), 2)) /
                                              (2.00 * Math.Pow(Funcs.v(prms, a[i]), 2)) +
                                              (2.00 * Math.Pow(Funcs.vp(prms, a[i]), 2) * Math.Pow(Funcs.mp(prms), 2)) /
                                              (Math.Pow(Funcs.v(prms, a[i]), 2) * Funcs.vpp(prms) -
                                               2.00 * Math.Pow(Funcs.vp(prms, a[i]), 2) * Funcs.v(prms, a[i])));
                                int[] h1Prtls = {-1, -1, i, -1};
                                double h1 = GetPNR.Run(type, prms, a, fxTD, rf0, n, nbuckets, prec, h1Prtls, plproc);
                                int[] h2Prtls = {-1, -1, -1, i};
                                double h2 = GetPNR.Run(type, prms, a, fxTD, rf0, n, nbuckets, prec, h2Prtls, plproc);
                                hess[i, j] = k1 * h1 + k2 * h2 + k3 * stdpnr;
                            }
                        }
                    }
                    else
                    {
                        hess[i, j] = hess[j, i];
                    }
                }
                if (i == fxTD - 1)
                {
                    Trace.Write(" (Done)");
                }
                Trace.WriteLine("");
            }
            Trace.WriteLine("");

            // Return the Hessian matrix.
            return hess;
        }
    }
}
