using System;
using System.Collections.Generic;
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
            if (args[0].EndsWith(".txt"))
            {
                Import(args[0]);
            }
            else
            {
                Extract(args[0]);
            }
            Console.WriteLine("Done!");
        }

        private static void Help()
        {
            string executableName = GetExecutableName();
            Console.WriteLine("Usage:");
            Console.WriteLine("  {0} <loc_file>", executableName);
            Console.WriteLine("  Export <loc_file> to <loc_file>.txt.");
            Console.WriteLine(" OR");
            Console.WriteLine("  {0} <loc_file>.txt", executableName);
            Console.WriteLine("  Import from <loc_file>.txt to <loc_file>. The LOC file must already exist.");
            Environment.Exit(42);
        }

        private static string GetExecutableName()
        {
            return Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location);
        }

        private static void Extract(string locFileName)
        {
            Loc loc = new Loc();
            using (MemoryMappedFile file = MemoryMappedFile.CreateFromFile(locFileName, FileMode.Open))
            using (Stream stream = file.CreateViewStream())
            using (MemoryMappedViewAccessor accessor = file.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read))
            {
                loc.Read(stream, accessor);
            }
            loc.DecodeStrings();
            loc.Export(locFileName);
        }

        private static void Import(string textFileName)
        {
            Dictionary<uint, string> newEntries = new Dictionary<uint, string>();
            using (StreamReader text = new StreamReader(textFileName, System.Text.Encoding.Unicode))
            {
                while (true)
                {
                    string line = text.ReadLine();
                    if (line == null)
                    {
                        break;
                    }
                    if (line.Trim() == "")
                    {
                        continue;
                    }
                    int equalsPos = line.IndexOf('=');
                    uint key = uint.Parse(line.Substring(0, equalsPos));
                    string value = line.Substring(equalsPos + 1).Replace("[CR]", "\r").Replace("[LF]", "\n");
                    newEntries[key] = value;
                }
            }
            string locFileName = textFileName.Substring(0, textFileName.Length - 3);
            Loc loc = new Loc();
            using (MemoryMappedFile file = MemoryMappedFile.CreateFromFile(locFileName, FileMode.Open))
            using (Stream stream = file.CreateViewStream())
            using (MemoryMappedViewAccessor accessor = file.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read))
            {
                loc.Read(stream, accessor);
            }
            loc.DecodeStrings();
            Console.WriteLine("Writing new LOC file...");
            loc.Update(newEntries);
            loc.ReCompress();
            using (Stream stream = new FileStream(locFileName, FileMode.Create))
            {
                loc.Write(stream);
            }
        }
    }
}
