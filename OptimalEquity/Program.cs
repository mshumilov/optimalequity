/*
/
/ Inline functions:
/
/ 1.) m() = The real/expense-adjusted return mean as a function of the equity ratio at a given time-point.
/ 2.) v() = The real/expense-adjusted return variance as a function of the equity ratio at a given time-point.
/ 3.) mp() = The derivative of m() with respect to the equity ratio at a given time-point.
/ 4.) vp() = The derivative of v() with respect to the equity ratio at a given time-point.
/ 5.) vpp() = The 2nd derivative of v() with respect to the equity ratio at a given time-point.
/ 6.) mva() = The equity ratio that reflects the minimum variance portfolio given the stock/bond return variance & covariance.
/ 7.) kh1() = A quantity used in defining the special Hessian density h1().
/ 8.) f() = The density that represents real/expense-adjusted returns (historical distribution).
/ 9.) g() = The special density used for computing the gradient at each time-point, and also used to construct the Hessian.
/ 10.) h1() = The 1st special density used for computing the diagonal Hessian elements at each time-point.
/ 11.) h2() = The 2nd special density used for computing the diagonal Hessian elements at each time-point.
/
/ Function prototypes:
/
/ 1.) GetPNR() = Prototype for the function that returns the probability of avoiding ruin (i.e., success probability).
/ 2.) GetConst() = Prototype for the function that returns the vector of constants used to build any CDF call within this application.
/ 3.) GetGamma() = Prototype for the function that returns the vector of gamma RVs used to build any CDF call within this application.
/ 4.) GetCDF() = Prototype for the function that returns any CDF value needed within this application using the vectors of constants and gamma RVs.
/ 5.) ThrdPNRsim() = Prototype for function that threads the simulation ruin computation across CPU cores to reduce processing time.
/ 6.) ThrdPNRdyn() = Prototype for function that threads the dynamic programming ruin computation across CPU cores to reduce processing time.
/ 7.) PNRsim() = Prototype for function that derives the probability of no ruin (using simulation) for a given glide-path/withdrawal rate
/ 8.) PNRdyn() = Prototype for function that derives the probability of no ruin (using dynamic programming) for a given glide-path/withdrawal rate.
/ 9.) BldGrad() = Prototype for function that accepts an empty array and populates it with the gradient vector elements for a given glide-path.
/ 10.) Climb() = Prototype for function that climbs in the direction of the gradient and stops when no further progress can be made.
/ 11.) DrvHess() = Prototype for function that derives and returns the Hessian matrix for a given glide-path.
/ 12.) WrtAry() = Prototype for function that writes an array (either the glide-path or the gradient) to standard output or a file as a block.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OptimalEquity.Helpers;

namespace OptimalEquity
{
    class Program
    {
        static void Main(string[] args)
        {
            AppHelper.EnableLogging();
            AppHelper.UseDotAsDecimalSeparatorInStrings();

            var startTime = DateTime.Now;

            try
            {
                if (args.Length > 1)
                {
                    Trace.WriteLine(
                        "ERROR: Parameter misspecification. Incorrect # of parameters to the executable (expecting zero or one...).");
                    Trace.WriteLine("EXITING...main()...");
                    Console.Read();
                    Environment.Exit(1);
                }

                var concurrency = args.Length == 1 ? int.Parse(args[0]) : 0;
                // When concurrency==0, replace with the # of independent processing units on the computer running the application.
                // Note that main() runs in its own thread but that is not accounted for here as it sits idle.
                if (concurrency == 0)
                    concurrency = Environment.ProcessorCount + 1;
                if (concurrency == 1)
                    concurrency++; // If just one processing unit use 2 threads.

                AppHelper.InitializeInputFilesIfNotExists();

                Model.StaticGP.Run(concurrency);

            }
            finally
            {
                Trace.WriteLine("");
                Trace.WriteLine($"Time spent: {(DateTime.Now - startTime).TotalSeconds:F1} seconds");
                Trace.WriteLine("");
                Trace.Write("Done (hit return to exit).");
                Console.Read();
            }
        }
    }
}
