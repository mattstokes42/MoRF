//Matthew E. Stokes
//This file contains the overall routine for running MoRF on a group of datasets

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

class MoRF_Runner
{
    public static void Main(string[] args)
    {
        //Parse command line arguments, if any
        CommandLineUtils.Arguments parameters = new CommandLineUtils.Arguments(args);
        string method = "SWRFstar";
        if (parameters["f"]!=null)
            method = parameters["f"];
        int samples = 200; int attributes = 1000;
        if (parameters["s"]!=null)
            samples = Convert.ToInt16(parameters["s"]);
        int startMod = 0; int endMod = 69;
        if (parameters["sm"] != null && parameters["em"] != null)
        {
            startMod = Convert.ToInt16(parameters["sm"]);
            endMod = Convert.ToInt16(parameters["em"]);
            if (startMod < 0 || startMod > 69)
            {
                Console.WriteLine("Invalid startMod, resetting to 0");
                startMod = 0;
            }
            if (endMod < startMod || endMod > 69)
            {
                Console.WriteLine("Invalid endMod, resetting to 69");
                endMod = 69;
            }
        }
        int startRep = 0; int endRep = 99;
        if (parameters["sr"] != null && parameters["er"] != null)
        {
            startRep = Convert.ToInt16(parameters["sr"]);
            endRep = Convert.ToInt16(parameters["er"]);
            if (startRep < 0 || startRep > 99)
            {
                Console.WriteLine("Invalid startRep, resetting to 0");
                startRep = 69;
            }
            if (endRep < startRep || endRep > 99)
            {
                Console.WriteLine("Invalid endRep, resetting to 99");
                endRep = 99;
            }
        }
        List<double> nwParameters = new List<double> { };
        if (parameters["p1"] != null)
            nwParameters.Add(Convert.ToDouble(parameters["p1"]));
        if (parameters["p2"] != null)
            nwParameters.Add(Convert.ToDouble(parameters["p2"]));

        Console.WriteLine("Using "+ method +" weighting for models " +startMod +" to " +endMod+
            " with reps "+ startRep + " to " + endRep);
        //Run Relief on up 70 models with up to 100 replications each
        //Command line flags -sm (startMod) -em (endMod) -sr (startRep) -er (endRep)
        for (int m = startMod; m <= endMod; m++)
        { 
            string model= m.ToString();
            if (m < 10)
              model = "0" + m.ToString();
            for (int rep = startRep; rep <= endRep; rep++)
            {
                string filename = samples+"//"+model+"//"+ model + "." + samples + "." + rep + ".txt";
                //Array R holds the data - rows are samples, columns are SNPs, last column is class
                if (System.IO.File.Exists(filename))
                {
                    int[][] R = FileUtilities.Utils.ReadDataFromFile(filename);
                    samples = R.GetLength(0); attributes = R[0].GetLength(0);
                    //dW contains prior knowledge attribute weighting
                    double[] dW = new double[attributes];
                    for (int i = 0; i < attributes; i++)
                        dW[i] = 1;               
                    //This line runs the Relief algorithm on data R
                    double[] W = ReliefUtils.ReliefUtils.MoRFScore(method, R, dW,nwParameters);
                    string outputFilename = method + ".";
                    for (int p = 0; p < nwParameters.Count; p++)
                        outputFilename = outputFilename + nwParameters[p] + ".";
                    outputFilename = outputFilename + samples + "." + model + "." + rep + ".txt";
                    FileUtilities.Utils.WriteReliefOutput(outputFilename, W);
                }
                else Console.WriteLine("File " + filename + " not found!");
            }
          //Once done with all the reps, create summary statistics (power)
            Console.WriteLine("Generating summary statistics for model "+model+"...");
            string filenameStub = method + ".";
                for (int p = 0; p < nwParameters.Count; p++)
                    filenameStub = filenameStub + nwParameters[p] + ".";
            filenameStub = filenameStub + samples + "." + model + ".";
            FileUtilities.Utils.powerOutput(filenameStub, m, startRep, endRep);
        }
    }

}