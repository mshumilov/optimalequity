using System;

namespace OptimalEquity.Helpers
{
    public static class Funcs
    {
        public static double m(double[] prms, double a)
        {
            return (1.00 - prms[5]) * (1 + a * prms[0] + (1.00 - a) * prms[2]);
        }
        public static double v(double[] prms, double a)
        {
            return Math.Pow(1.00 - prms[5], 2) * (a * a * prms[1] + (1.00 - a) * (1.00 - a) * prms[3] + 2.00 * a * (1.00 - a) * prms[4]);
        }
        public static double mp(double[] prms)
        {
            return (1.00 - prms[5]) * (prms[0] - prms[2]);
        }
        public static double vp(double[] prms, double a)
        {
            return Math.Pow(1.00 - prms[5], 2) * (2.00 * a * prms[1] - 2.00 * (1.00 - a) * prms[3] + (2.00 - 4.00 * a) * prms[4]);
        }
        public static double vpp(double[] prms)
        {
            return Math.Pow(1.00 - prms[5], 2) * (2.00 * prms[1] + 2.00 * prms[3] - 4.00 * prms[4]);
        }
        public static double mva(double[] prms)
        {
            return (prms[3] - prms[4]) / (prms[1] + prms[3] - 2.00 * prms[4]);
        }
        public static double kh1(double[] prms, double a)
        {
            return (-2.00 * vp(prms, a) * mp(prms) * v(prms, a)) / (v(prms, a) * vpp(prms) - 2.00 * Math.Pow(vp(prms, a), 2));
        }
        public static double f(double[] prms, double a, double r)
        {
            return (1 / Math.Sqrt(2.00 * Math.PI * v(prms, a))) * Math.Exp(-Math.Pow(r - m(prms, a), 2) / (2.00 * v(prms, a)));
        }
        public static double g(double[] prms, double a, double r)
        {
            return f(prms, a, r) * Math.Pow(vp(prms, a) * (r - m(prms, a)) + mp(prms) * v(prms, a), 2) / (Math.Pow(mp(prms), 2) * Math.Pow(v(prms, a), 2) + v(prms, a) * Math.Pow(vp(prms, a), 2));
        }
        public static double h1(double[] prms, double a, double r)
        {
            return f(prms, a, r) * Math.Pow((r - m(prms, a) + kh1(prms, a)), 2) / (v(prms, a) + Math.Pow(kh1(prms, a), 2));
        }
        public static double h2(double[] prms, double a, double r)
        {
            return f(prms, a, r) * 2.00 * Math.Pow(Math.Pow(r - m(prms, a), 2) * vp(prms, a) / (2.00 * v(prms, a)) + mp(prms) * (r - m(prms, a)) - vp(prms, a) / 2.00, 2) / (Math.Pow(vp(prms, a), 2) + 2.00 * v(prms, a) * Math.Pow(mp(prms), 2));
        }
    }
}
