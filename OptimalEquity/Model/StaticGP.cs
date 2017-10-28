/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: StaticGP.cpp (Defines the entry point for the console application.)
/
/ Function: main()
/
/ Summary:
/
/ After launching the executable the user is prompted for the location of a directory where the control file and initial glide-path is stored in files
/ named "control.txt" and "gp.txt". The file names can be changed in the header. The control file contains the means/variances/covariance for stock and
/ bond returns (real) which are assumed normally distributed and constant with respect to time for this application. (Neither is a requirement for the
/ proposed method, only for this particular implementation.) The expense ratio, retirement horizon length, fixed inflation-adjusted withdrawal rate,
/ optimization method (gradient ascent/Newton-Raphson), convergence threshold, estimation method (dynamic program/simulation), base simulation sample
/ size and testing alphas, dynamic program discretization size and maximum ruin factor are also included in the control file. The glide-path file holds
/ the initial glide-path in a single column which is the starting point for the optimization. When using gradient ascent, the direction of steepest
/ ascent is computed for the current glide-path and climbing begins in that direction until no further progress can be made. When estimating gradient
/ elements using simulation, no progress is made when the null hypothesis that the new probability of success is >= the old probability of success is
/ rejected (using the 1st alpha from the control file). (Note that the probabilities being compared are subject to sampling error when estimating
/ probabilities using simulation which is why a hypothesis test is used to determine when to stop climbing.) When using dynamic programming to estimate
/ success probabilities we use a direct comparison, not a hypothesis test. When climbing progress ends we recompute the new gradient and proceed to
/ climb in the new direction until no further progress is made. This procedure is repeated until the largest absolute effective value of the gradient
/ vector is below the convergence threshold specified in the control file. By driving the gradient vector to zero in all dimensions we move to a point
/ on the surface that reflects a local/global maximum, which is true when the region is concave (i.e., Hessian matrix is negative semi-definite <=> all
/ its eigenvalues are non-positive.) If all initial glide-paths exist in a concave region and lead to the same local optimum we consider it an empirical
/ interior point optimum and argue against the need for a meta-heuristic. All elements of the gradient vector and Hessian matrix can be expressed as
/ probabilities (or differences/linear combinations of) and these are estimated using either dynamic programming or simulation. When climbing using
/ simulation a sample size of 2X the base sample size is used since probabilities are compared. Once no further progress is made we recompute the
/ success probability for that glide-path to remove any built-in upward sampling bias using a simulation sample size of 4X the base sample size.
/ Gradient and Hessian computations use the base sample size. This probability is used when computing the subsequent gradient. When using gradient
/ ascent with simulation we test each element of the gradient for equality with zero and set it to zero if the hypothesis test is not rejected (using the
/ 2nd alpha from the control file). Since these quantities are subject to sampling error we do not want to climb in the direction of that sampling error
/ if it can be avoided. The resulting vector is referred to as the adjusted gradient and only applies when simulation is the estimation method. The
/ Newton-Raphson optimization method is also available and uses a first-order Taylor expansion to estimate the gradient in the neighborhood of a given
/ point. We know that the optimal value is the point where the gradient equals zero and this approximation is set to zero and solved yielding a new
/ point. The new point becomes the old point and the process is repeated iteratively until the convergence criteria is met. Convergence is achieved for
/ both optimization methods when the largest absolute gradient entry is smaller than the value supplied in the control file. The user thus has 2 options
/ for optimization (gradient ascent or Newton-Raphson) and 2 options for estimating success probabilities (simulation or dynamic programming). The final
/ glide-path is written to the file output.txt located in the directory specified by the user. The Hessian and its eigenvalues are calculated at the
/ final glidepath to confirm that the region is concave.
/
/ Parameters:
/
/ The executable takes either one or no arguments when it is launched. This value reflects the number of concurrent processes to utilize when running
/ the program. If the user does not supply a value then the program will determine the maximum number allowed on the PC running it and use this value.
/
/ Input Files:
/
/ 1.) Control file named "control.txt" in the directory that the user provides. Change the name in the header file if desired. This value is a constant
/ set in the header to the constant paramfile.
/ 2.) Initial glide-path which is in the file specified by the header constant initgpfile (default "gp.txt"). This file should contain exactly TD values
/ between MVA and 1.00 where MVA is the minimum variance alpha for this arrangement. We do not use alphas below the MVA. This file is also located
/ in the directory specified by the user.
/
/ Output Files:
/
/ 1.) The file "output.txt" is written to the directory specified by the user and contains the final optimal glide-path, assuming convergence is
/ achieved.
/
/ Return Value:
/
/ A value of 0 (for success) or an appropriate error code (for failure).
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Accord.Math.Decompositions;

namespace OptimalEquity.Model
{
    public static class StaticGP
    {
        public static void Run(int pllproc)
        {
            double[] @params = new double[6];
            double maxgrad = 0;
            double stpthrsh = 0;
            double maxeigen;
            double mineigen;
            double wr = 0;
            double[] alpha = { 0.50, 1.00 };
            double rfMax = 2.75;
            int td = 0;
            int[] prtls = { -1, -1, -1, -1 };
            int iterint = 1;
            long sn = (long)Math.Pow(10.0, 7);
            int prec = (int)Math.Pow(10.0, 4);
            string type = "";
            string alg = "";

            const string rootdir = Config.ConfigsRootDir;

            // Read setup details from control file.
            try
            {
                var getParams = File.ReadAllText(Path.Combine(rootdir, Config.ControlFile));
                var s = getParams.Split(new[] { ' ', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

                @params = s.Take(6).Select(i => double.Parse(i, CultureInfo.InvariantCulture)).ToArray();
                td = int.Parse(s[6]);
                wr = double.Parse(s[7]);
                stpthrsh = double.Parse(s[8]);
                alg = s[9];
                type = s[10];
                if (type.Equals("sim"))
                {
                    sn = int.Parse(s[11]);
                    alpha[0] = double.Parse(s[12]);
                    alpha[1] = double.Parse(s[13]);
                }
                else if (type.Equals("dp"))
                {
                    prec = int.Parse(s[11]);
                    rfMax = double.Parse(s[12]);
                }
            }
            catch (Exception ex)
            {
                Trace.Write("ERROR: Could not read file: ");
                Trace.WriteLine(Path.Combine(rootdir, Config.ControlFile) + $". {ex.Message}");
                Trace.WriteLine("EXITING...main()...");
                Console.Read();
                Environment.Exit(1);
            }


            // Read initial glide-path from file.
            var gp = new double[td];

            var gpPath = Path.Combine(rootdir, Config.InitGladepathFile);
            try
            {
                gp = File.ReadAllLines(gpPath).Select(double.Parse)
                    .Take(td).ToArray();
            }
            catch (Exception ex)
            {
                Trace.Write("ERROR: Could not read file: ");
                Trace.WriteLine(gpPath + ". Error: " + ex.Message);
                Trace.WriteLine("EXITING...main()...");
                Console.Read();
                Environment.Exit(1);

            }

            if (gp.Length != td)
            {
                Trace.Write("ERROR: File: ");
                Trace.Write(gpPath);
                Trace.Write(" needs ");
                Trace.Write(td);
                Trace.WriteLine(" initial asset allocations, but has fewer.");
                Trace.WriteLine("EXITING...main()...");
                Console.Read();
                Environment.Exit(1);
            }

            // Display optimization algorithm.
            Trace.Write(@"===> Optimization algorithm: ");
            if (alg == "nr")
                Trace.WriteLine(@"Newton's Method");
            else if (alg == "ga")
                Trace.WriteLine(@"Gradient Ascent");

            // Display estimation method.
            Trace.WriteLine("");
            Trace.Write(@"===> Estimation method: ");
            if (type == "sim")
            {
                Trace.WriteLine(@"Simulation");
            }
            else if (type == "dp")
            {
                Trace.WriteLine(@"Dynamic Program");
                pllproc = 4 * pllproc;
            }

            // Declare variables that depend on data read from the control file for sizing.
            double[,] hess;
            double[] ngrdnt = new double[td];
            EigenvalueDecomposition hevals;

            var grad = new double[td];

            // Take some steps (1 full iteration but no more than 50 steps) in the direction of steepest ascent. This can move us off
            // the boundary region where computations may be unstable (infinite), especially when constructing the Hessian for Newton's method.
            // Also, this initial stepping usually makes improvements very quickly before proceeding with the optimization routine.
            double probnr = GetPNR.Run(type, @params, gp, td, wr, 4 * sn, (int)(rfMax * prec), prec, prtls, pllproc);
            Trace.WriteLine("");
            Trace.WriteLine("Initial Glide-Path (w/Success Probability):");
            WrtAry.Run(probnr, gp, "GP", td);
            for (int s = 1; s <= 2; ++s)
            {
                maxgrad = BldGrad.Run(type, @params, gp, td, wr, sn, (int)(rfMax * prec), prec, probnr, 4 * sn, alpha[1], pllproc, grad);
                if (maxgrad <= stpthrsh)
                {
                    Trace.Write("The glide-path supplied satisfies the EPSILON convergence criteria: ");
                    Trace.WriteLine($"{maxgrad:F15} vs. {stpthrsh:F15}");
                    s = s + 1;
                }
                else if (s != 2)
                {
                    probnr = Climb.Run(type, @params, gp, td, wr, 2 * sn, (int) (rfMax * prec), prec, pllproc, maxgrad,
                        probnr, 4 * sn, grad, alpha[0], 50);
                    Trace.WriteLine("");
                    Trace.WriteLine("New (Post Initial Climb) Glide-Path (w/Success Probability):");
                    WrtAry.Run(probnr, gp, "GP", td);
                }
                else if (maxgrad <= stpthrsh)
                {
                    Trace.Write("The glide-path supplied satisfies the EPSILON convergence criteria after intial climb without iterating: ");
                    Trace.WriteLine($"{maxgrad:F15} vs. {stpthrsh:F15}");
                }
            }

            // Negate the gradient if using NR method.
            if (alg == "nr")
            {
                for (int y = 0; y < td; ++y)
                    ngrdnt[y] = -1.00 * grad[y];
            }

            // If convergence is not achieved after initial climb then launch into full iteration mode.
            while (maxgrad > stpthrsh)
            {
                Trace.WriteLine("");
                Trace.WriteLine("=========================");
                Trace.WriteLine($"Start Iteration #{iterint}");
                Trace.WriteLine("=========================");
                if (alg == "nr")
                {
                    // Record the probability before iterating.
                    double strtpnr = probnr;

                    // Build the Hessian matrix for this glide-path and derive its eigenvalues. (Display the largest & smallest value.)
                    // This is required when method=nr. When either procedure ends with convergence we recompute the Hessian matrix to
                    // ensure we are at a local/global maximum (done below after convergence).
                    hess = DrvHess.Run(type, @params, gp, td, wr, sn, (int)(rfMax * prec), prec, pllproc, grad, probnr);
                    //hevals.compute(hess, false);
                    hevals = new EigenvalueDecomposition(hess);

                    var reals = hevals.RealEigenvalues;
                    maxeigen = reals.Max();
                    mineigen = reals.Min();

                    // Display the smallest/largest eigenvalues.
                    Trace.WriteLine("");
                    Trace.Write("Min Hessian eigenvalue for this iteration (>=0.00 --> convex region): ");
                    Trace.WriteLine(mineigen);
                    Trace.WriteLine("");
                    Trace.Write("Max Hessian eigenvalue for this iteration (<=0.00 --> concave region): ");
                    Trace.WriteLine(maxeigen);

                    // Update the glidepath and recompute the probability using the new glidepath.
                    //sol = hess.colPivHouseholderQr().solve(ngrdnt);
                    var qr = new QrDecomposition(hess);
                    var sol = qr.Solve(ngrdnt);

                    for (int y = 0; y < td; ++y)
                    {
                        gp[y] += sol[y];
                    }
                    probnr = GetPNR.Run(type, @params, gp, td, wr, 4 * sn, (int)(rfMax * prec), prec, prtls, pllproc);

                    // If success probability has worsened alert the user.
                    if (probnr < strtpnr)
                    {
                        Trace.WriteLine("");
                        Trace.WriteLine("NOTE: The success probability has worsened during the last iteration. This could happen for different reasons:");
                        Trace.WriteLine(" 1.) The difference in probabilities is beyond the system's ability to measure accurately (i.e., beyond 15 significant digits).");
                        Trace.WriteLine(" 2.) The difference is due to estimation/approximation error.");
                        Trace.WriteLine(" 3.) You may be operating along the boundary region. In general the procedure is not well defined on the boundaries. (Try gradient ascent.)");
                    }
                }
                else if (alg == "ga")
                {
                    // Update the glide-path and recompute the probability using the new glide-path.
                    probnr = Climb.Run(type, @params, gp, td, wr, 2 * sn, (int) (rfMax * prec), prec, pllproc, maxgrad,
                        probnr, 4 * sn, grad, alpha[0]);
                }

                // Display the new glide-path.
                Trace.WriteLine("");
                Trace.Write("New Glide-Path:");
                WrtAry.Run(probnr, gp, "GP", td);

                // Rebuild the gradient and negate it when using NR.
                maxgrad = BldGrad.Run(type, @params, gp, td, wr, 1 * sn, (int) (rfMax * prec), prec, probnr,
                    4 * sn, alpha[1], pllproc, grad);
                if (alg == "nr")
                {
                    for (int y = 0; y < td; ++y)
                        ngrdnt[y] = -1.00 * grad[y];
                }
                // Report the convergence status.
                Trace.WriteLine("");
                Trace.WriteLine($"EPSILON Convergence Criteria: {maxgrad:F15} vs. {stpthrsh:F15}");
                if (maxgrad <= stpthrsh)
                {
                    Trace.WriteLine("");
                    Trace.WriteLine("==========> EPSILON Convergence criteria satisfied. <==========");
                }
                Trace.WriteLine("");
                Trace.WriteLine(new String('=', 25));
                Trace.Write("End Iteration #");
                Trace.WriteLine(iterint);
                Trace.WriteLine(new String('=', 25));
                iterint++;
            }

            // Build Hessian and confirm we are at a maximum, not a saddle-point or plateau for example.
            Trace.WriteLine("");
            Trace.WriteLine("Convergence Achieved: Final step is to confirm we are at a local/global maximum. Hessian is being built.");
            hess = DrvHess.Run(type, @params, gp, td, wr, sn, (int)(rfMax * prec), prec, pllproc, grad, probnr);

            hevals = new EigenvalueDecomposition(hess);
            var r = hevals.RealEigenvalues;
            maxeigen = r.Max();
            mineigen = r.Min();

            // Display the smallest/largest eigenvalues.
            Trace.WriteLine("");
            Trace.Write("Min Hessian eigenvalue at solution [>=0.00 --> convex region --> (local/global) minimum]: ");
            Trace.WriteLine(mineigen);
            Trace.WriteLine("");
            Trace.Write("Max Hessian eigenvalue at solution [<=0.00 --> concave region --> (local/global) maximum]: ");
            Trace.WriteLine(maxeigen);

            // Write final GP to the output file.
            Trace.WriteLine("");
            if (maxeigen <= 0 || mineigen >= 0)
                Trace.Write("(Local/Global) Optimal ");
            Trace.WriteLine("Glide-Path:");
            WrtAry.Run(probnr, gp, "GP", td, Path.Combine(rootdir, Config.Outfile));
            Trace.WriteLine("");
        }
    }
}
