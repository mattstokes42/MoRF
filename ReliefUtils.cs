//Matthew E. Stokes
//This file contains routines for computing the MoRF score of each attribute in a given dataset

using System;
using System.Collections.Generic;
using System.Text;

namespace ReliefUtils
{
    public class ReliefUtils
    {
        //This function returns the "difference value" between the same attribute on two individuals
        //Function "diff" in paper nomenclature
        public static double attributeDifference(int[][] R, int attribute, int inst1, int inst2)
        {
            if (R[inst1][attribute] == R[inst2][attribute])
                return 0;
            else return 1;
        }
        //This function returns a "difference score" matrix of pairwise distances between
        //the whole sequence of two individuals (high scores means less similarity)
        public static double[,] distanceMeasure(int[][] R, double[] dW)
        {
            int attributes = R[0].GetLength(0); int samples = R.GetLength(0);
            double[,] pairwiseDistance = new double[samples,samples];
            for (int i = 0; i < samples; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    //This is the weighted edit distance
                    for (int a = 0; a < attributes; a++)
                        pairwiseDistance[i, j] += dW[a] * attributeDifference(R, a, i, j);
                }
            }
            return pairwiseDistance;
        }

        //Function "c" in the paper nomenclature
        public static double classComparator(int class1, int class2)
        {
            if (class1 == class2)
                return -1;
            else return 1;
        }

        //Find the mean pairwise distance to pick a neighbor weighting function center
        public static double calculateT(double[,] score)
        {
            int samples = score.GetLength(0);
            double T = 0.0; int n = 0;
            for (int i = 0; i < samples; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    if (i != j)
                    {
                        T += score[i, j];
                        n++;
                    }
                }
            }
            return T / n;
        }

        //Find the standard deviation of pairwise distances
        public static double calculateS(double[,] score, double mean)
        {
            int samples = score.GetLength(0);
            double S = 0.0; int n = 0;
            for (int i = 0; i < samples; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    if (i != j)
                    {
                        S += Math.Pow(score[i, j] - mean, 2.0);
                        n++;
                    }
                }
            }
            return Math.Sqrt(S / (n + 1));
        }

        //These are the neighbor weighting functions (command-line flag "-f")
        //Function "f" in paper nomenclature
        public static double[,] neighborWeight(string method, double[,] distanceMat, int[][] R, List<double> parameter)
        {
            int samples = distanceMat.GetLength(0);
            double[,] f = new double [samples,samples];
            double T = calculateT(distanceMat);
            double S = calculateS(distanceMat, T);
            for (int i = 0; i < samples; i++)
            {
                //ReliefF weighting finds n nearest hits and misses to each query sample
                //parameter[0] = k, the number of nearest neighbors to use
                if (method == "ReliefF" || method == "ReliefFstar")
                {
                    List<sample> hits = new List<sample> { };
                    List<sample> misses = new List<sample> { };
                    //Put all neighbor samples in lists of Hits and Misses
                    for (int j = 0; j < samples; j++)
                    {
                        if (j != i)
                        {
                            if (R[i][R[0].GetLength(0) - 1] == R[j][R[0].GetLength(0) - 1])
                                hits.Add(new sample(j, distanceMat[i, j]));
                            else
                                misses.Add(new sample(j, distanceMat[i, j]));
                        }
                    }
                    //Sort hits and misses by ascending distance from query sample
                    hits.Sort(delegate(sample s1, sample s2) { return s1.score.CompareTo(s2.score); });
                    misses.Sort(delegate(sample s1, sample s2) { return s1.score.CompareTo(s2.score); });
                    //Take k nearest hits and k nearest misses, ignore the rest
                    int k = 10;
                    if (parameter.Count!=0)
                        k = Convert.ToInt16(parameter[0]);
                    for (int t = 0; t < k; t++)
                    {
                        f[i, hits[t].ID] = 1;
                        f[i, misses[t].ID] = 1;
                    }
                    //ReliefFstar gives opposite weight to k farthest hits and misses as well
                    if (method == "ReliefFstar")
                    {
                        for (int n = 0; n < k; n++)
                        {
                            f[i, hits[hits.Count-1-n].ID] = -1;
                            f[i, misses[misses.Count-1-n].ID] = -1;
                        }
                    }
                }
                else
                {
                    //Spatially sensitive versions use global, not sample-wise statistics
                    for (int j = 0; j < samples; j++)
                    {
                        if (method == "SWRFstar")
                        {
                            //Sigmoid weighting function
                            //parameter[0] is the scaling factor of SWRFstar
                            double scale = 4.0;
                            if (parameter.Count != 0)
                                scale = parameter[0];
                            f[i, j] = 2 * (1 / (1 + Math.Pow(Math.E, -(T - distanceMat[i, j]) / (S / scale)))) - 1;
                        }
                        if (method == "SWRF")
                        {
                            //Half-sigmoid weighting function
                            //parameter[0] is the scaling factor of SWRF
                            double scale = 4.0;
                            if (parameter.Count != 0)
                                scale = parameter[0];
                            if (distanceMat[i, j] < T)
                                f[i, j] = 2 * (1 / (1 + Math.Pow(Math.E, -(T - distanceMat[i, j]) / (S / scale)))) - 1;
                        } 
                        if (method == "SURF")
                        {
                            //Step weighting function
                            if (distanceMat[i, j] < T)
                                f[i, j] = 1;
                            else
                                f[i, j] = 0;
                        }
                        if (method == "SURFstar")
                        {
                            //SURF* weighting function
                            if (distanceMat[i, j] < T)
                                f[i, j] = 1;
                            else
                                f[i, j] = -1;
                        }
                        if (method == "SURF-S")
                        {
                            //Nearer than T-S step weighting function
                            if (distanceMat[i, j] < (T - S))
                                f[i, j] = 1;
                        }
                        if (method == "ramp")
                        {
                            if (distanceMat[i,j]<T)
                                f[i, j] = 1 - (distanceMat[i, j]/T);
                        }
                        if (method == "rampstar")
                        {
                                f[i, j] = 1 - (distanceMat[i, j] / T);
                        }
                    }
                }           
                f[i, i] = 0;
            }
            return f;
        }

        public static double[] MoRFScore(string method, int[][] R, double[] dW,List<double> parameters)
        {
            int samples = R.GetLength(0);
            int attributes = R[0].GetLength(0)-1;
            double[] W = new double[attributes];
            //Calculate pairwise distance matrix and neighbor weight matrix
            double[,] pairwiseDist = distanceMeasure(R, dW);
            double[,] nW = neighborWeight(method, pairwiseDist, R,parameters);
            Console.WriteLine("Calculating MoRF score using "+method+" neighbor weighting...");
            //Calculate scores for all attributes simultaneously by pairwise sample comparisons
            //Each query sample has a contribution (delta W) to the overall score of each attribute
            for (int querySample = 0; querySample < samples; querySample++)
            {
                double[] num = new double[attributes]; 
                double[] den = new double[attributes];
                //Examine each neighbor of the query sample
                for (int i = 0; i < samples; i++)
                {
                    if (i != querySample)
                    {
                        double f = nW[querySample,i];
                        double c = classComparator(R[querySample][attributes],R[i][attributes]);
                        for (int a = 0; a < attributes; a++)
                        {
                            double diff = attributeDifference(R, a, querySample, i);
                            num[a] += c * diff * f * dW[a];
                            den[a] += Math.Abs(f * dW[a]);
                        }
                    }
                }
                //Find the query sample's total contribution to the MoRF score of all atributes
                for (int a = 0; a < attributes; a++)
                    W[a] += num[a] / (den[a] * samples);
            }
            Console.WriteLine("Finished calcualting scores");
            return W;
        }
    }

    //Simple class to facilitate sorting SNP IDs by MoRF score
    public class SNP
    {
        public int ID;
        public double W;
        public SNP(int number, double weight)
        {
            ID = number;
            W = weight;
        }   
    }

    //Simple class to facilitate sorting neighbor samples by distance from query sample
    public class sample
    {
        public int ID;
        public double score;
        public sample(int number, double s)
        {
            ID = number;
            score = s;
        }
    }
}
