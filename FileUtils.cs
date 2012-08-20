//Matthew E. Stokes
//This file contains utilities for reading and writing data files

using System;
using System.IO;
using System.Collections.Generic;
using System.Text;

namespace FileUtilities
{
    public class Utils
    {
        //Reads data from tab-delimited file with header, last column is class
        public static int[][] ReadDataFromFile(string filename)
        {
            StreamReader SR = File.OpenText(filename);
            char delimiter = '\t';
            //Count entries in the header line to get # of SNPs
            int attributes = SR.ReadLine().Split(delimiter).Length-1;
            //Count lines to get the # of samples
            int samples = 0;
            while (SR.ReadLine()!=null)
                samples++;
            Console.WriteLine("Reading "+samples+" samples with "+attributes+" SNPs from "+filename+"...");
            SR.Close();
            //Now make the array and read in data
            int[][] R = new int[samples][];
            for (int i=0;i<samples;i++)
                R[i] = new int[attributes+1];
            SR = File.OpenText(filename);
            //Read off header line
            SR.ReadLine();
            int sampleNumber = 0;
            while (sampleNumber < samples)
            {
                //Read SNP values
                string[] parts = SR.ReadLine().Split(delimiter);
                for (int i = 0; i < R[sampleNumber].Length - 1; i++)
                    R[sampleNumber][i] = Convert.ToInt32(parts[i]);
                //Read class value
                R[sampleNumber][R[sampleNumber].Length - 1] = Convert.ToInt32(parts[parts.Length - 1]);
                sampleNumber++;
            }
            SR.Close();
            Console.WriteLine("Finished reading...");
            return R;
        }

        //Write the MoRF scores with a SNP ID and score on each line, separated by a tab
        public static void WriteReliefOutput(string filename, double[] W)
        {
            StreamWriter output = new System.IO.StreamWriter(filename);
            for (int i = 0; i < W.Length; i++)
                output.WriteLine(i+"\t"+W[i]);
            output.Close();
        }

        //Read the MoRF score output files (for power statistics)
        public static List<ReliefUtils.SNP> ReadReliefOutput(string filename)
        {
            List<ReliefUtils.SNP> SNPs = new List<ReliefUtils.SNP> { };
            StreamReader SR = File.OpenText(filename);
            string line = SR.ReadLine();
            while (line != null)
            {
                string[] parts = line.Split('\t');
                SNPs.Add(new ReliefUtils.SNP(Convert.ToInt16(parts[0]),Convert.ToDouble(parts[1])));
                line = SR.ReadLine();
            }
            return SNPs;
        }

        //Write the power output file
        //Output is 100 numbers separated by tabs, representing the power at 1% threshold, 2%, 3%... 100%
        public static void powerOutput(string filenameStub,int m, int startRep, int endRep)
        {
            int[] ranks = new int[endRep - startRep + 1];
            double[] power = new double[100];
            int attributes = 1000;
            for (int rep = startRep; rep <= endRep; rep++)
            {
                string filename = filenameStub + rep + ".txt";
                List<ReliefUtils.SNP> SNPs = FileUtilities.Utils.ReadReliefOutput(filename);
                attributes = SNPs.Count;
                //Sort SNPs by descending score
                SNPs.Sort(delegate(ReliefUtils.SNP S1, ReliefUtils.SNP S2) { return S2.W.CompareTo(S1.W); });
                //Find the lower ranked of the 2 informative SNPs (0 and 1)
                int rank0 = SNPs.FindIndex(delegate(ReliefUtils.SNP s) { return s.ID == 0; });
                int rank1 = SNPs.FindIndex(delegate(ReliefUtils.SNP s) { return s.ID == 1; });
                ranks[rep-startRep] = Math.Max(rank0, rank1);
            }
            string powerFilename = filenameStub + startRep + "to" + endRep + ".power.txt";
            StreamWriter powerOutput = new StreamWriter(powerFilename);
            //Find the power by percentiles (% of reps with both SNPs above threshold)            
            for (int percentile = 1; percentile <= 100; percentile++)
            {
                int passed = 0;
                for (int rep = 0; rep <= endRep - startRep; rep++)
                {
                    if (ranks[rep] <= (attributes * percentile / 100.0))
                        passed++;
                }
                power[percentile - 1] = passed / (endRep - startRep + 1.0);
                Console.WriteLine("Power at " + percentile + "%: " + power[percentile - 1]);
                powerOutput.Write(power[percentile - 1]);
                if (percentile != 100) powerOutput.Write("\t");
            }
            powerOutput.Close();
        }
    }
}
