using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {
            string fileN = "";
            while (true)
            {
               fileN= Console.ReadLine();
                int counter = 0;
                string line;
                string otto = "";

                // Read the file and display it line by line.
                System.IO.StreamReader file =
                    new System.IO.StreamReader(@"c:\"+fileN);
                while ((line = file.ReadLine()) != null)
                {
                    otto = otto + line;
                   // System.Console.WriteLine(line);
                    counter++;
                }

                file.Close();
          //      System.Console.WriteLine("There were {0} lines.", counter);
                int count = 1;
                foreach (char c in otto)
                {
                    if (c == ',')
                    {
                        count++;
                    }
                }
                Console.Write(fileN+" "+count+"\n");
            }
// Suspend the screen.
System.Console.ReadLine();
        }
    }
}
