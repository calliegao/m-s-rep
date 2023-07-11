# m-s-rep
Multiscale skeletal representation

This repository storges code and results for the paper "MR-based spatiotemporal anisotropic atrophy evaluation of hippocampus in AD progression by multiscale skeletal representation". 


Files in this repository:


SpokeInterpolation.py: Interpolate original s-rep.

RefineSpokes.py: Refine the spokes to fit precisely to the boundary surface of hippocampus.

forms: Morphological measurements extracted from bilateral hippocampi. 

StatMeasurementsGLM.R: Test the measurements using GLM. The input data is in the form file.

shapeStat.rar: Test the global shape difference represented by s-rep data.

     shapeStat\Vtp2SReppData.m: Normalize the s-rep data and arrange the vtp files of original s-rep data by columns.
     
     shapeStat\CPNS\my_calculateCPNS.m: Calculate CPNS mean s-rep of a group of subjects.
     
     shapeStat\CPNS\CalculatePermutationMeans.m: Step 1 of global shape testing, calculating means of each permutations. 
     
     shapeStat\CPNS\PermutationTest.m: Step 2 of globalshape testing, calculating statistics and p values.
