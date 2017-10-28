using System;
using System.Diagnostics;
using System.IO;
using OptimalEquity.Model;
using OptimalEquity.Properties;

namespace OptimalEquity.Helpers
{
    public static class AppHelper
    {
        /// <summary>
        /// Should be called once on app start.
        /// </summary>
        public static void EnableLogging()
        {
            if (!Directory.Exists(Config.ConfigsRootDir))
                Directory.CreateDirectory(Config.ConfigsRootDir);

            Trace.Listeners.Clear();
            TextWriterTraceListener twtl =
                new TextWriterTraceListener(Path.Combine(Config.ConfigsRootDir, AppDomain.CurrentDomain.FriendlyName + ".log"))
                {
                    Name = "TextLogger",
                    TraceOutputOptions = TraceOptions.ThreadId | TraceOptions.DateTime
                };
            ConsoleTraceListener ctl = new ConsoleTraceListener(false) {TraceOutputOptions = TraceOptions.DateTime};

            Trace.Listeners.Add(twtl);
            Trace.Listeners.Add(ctl);
            Trace.AutoFlush = true;
        }


        public static void UseDotAsDecimalSeparatorInStrings()
        {
            System.Globalization.CultureInfo customCulture =
                (System.Globalization.CultureInfo) System.Threading.Thread.CurrentThread.CurrentCulture.Clone();
            customCulture.NumberFormat.NumberDecimalSeparator = ".";
            System.Threading.Thread.CurrentThread.CurrentCulture = customCulture;
        }

        /// <summary>
        /// Create default input files if they are not created yet.
        /// </summary>
        public static void InitializeInputFilesIfNotExists()
        {
            if (!Directory.Exists(Config.ConfigsRootDir))
                Directory.CreateDirectory(Config.ConfigsRootDir);

            InitializeGlidepathFileIfNotExists();
            InitializeControlFileIfNotExists();
        }

        private static void InitializeGlidepathFileIfNotExists()
        {
            var filePath = Path.Combine(Config.ConfigsRootDir, Config.InitGladepathFile);
            if (File.Exists(filePath))
                return;
            File.WriteAllText(filePath, Resource.gp);
        }

        private static void InitializeControlFileIfNotExists()
        {
            var filePath = Path.Combine(Config.ConfigsRootDir, Config.ControlFile);
            if (File.Exists(filePath))
                return;
            File.WriteAllText(filePath, Resource.control);
        }
    }
}
