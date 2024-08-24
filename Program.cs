using System;
using System.IO;
using System.IO.MemoryMappedFiles;

namespace watch_dogs_loc
{
    class MainClass
    {
        public static void Main(string[] args)
        {
            Console.WriteLine("Watch_Dogs 1, 2, and Legion .loc tool - celikeins - 2019.06.07");
            if (args.Length != 1)
            {
                Help();  
            }
            Loc loc = new Loc();
            using (MemoryMappedFile file = MemoryMappedFile.CreateFromFile(args[0], FileMode.Open))
            using (Stream stream = file.CreateViewStream())
            using (MemoryMappedViewAccessor accessor = file.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read))
            {
                loc.Read(stream, accessor);
            }
            loc.Export(args[0]);
            Console.WriteLine("Done!");
        }

        private static void Help()
        {
            Console.WriteLine("Usage:");
            Console.WriteLine("  {0} <loc_file>", GetExecutableName());
            Console.WriteLine("  Export <loc_file> to <loc_file>.txt.");
            Environment.Exit(42);
        }

        private static string GetExecutableName()
        {
            return Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location);
        }
    }
}
