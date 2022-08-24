# Iterative Poisson Surface Reconstruction (iPSR) for Unoriented Points

Iterative Poisson Surface Reconstruction (iPSR) for Unoriented Points  
[Fei Hou](https://lcs.ios.ac.cn/~houf/),    Chiyu Wang,    [Wencheng Wang](https://lcs.ios.ac.cn/~whn/),    [Hong Qin](https://www3.cs.stonybrook.edu/~qin/),    Chen Qian,    [Ying He](https://personal.ntu.edu.sg/yhe/)  
*ACM Transactions on Graphics, 41, 4, Article 128, 13 pages, 2022. (SIGGRAPH 2022)*

iPSR extends the popular Poisson Surface Reconstruction ([https://github.com/mkazhdan/PoissonRecon](https://github.com/mkazhdan/PoissonRecon)). iPSR has no more need of oriented normals as input, but infers the normals in an iterative manner. It is used to reconstruct surface from only points input.  
[project page](https://lcs.ios.ac.cn/~houf/pages/ipsr/index.html)

### Compilation:
Windows: The code is tested by Visual Studio. The ipsr.vcxproj is an example to configure the project.  
Linux: The code is tested by GCC and Clang with makefile.  
Executable: [Win64](https://lcs.ios.ac.cn/~houf/pages/ipsr/iPSR.zip)

### Usage:
\-\-in &lt;input ply file name&gt;  
The input file should be in .ply format and only 3D point coordinates are needed.

\-\-out &lt;output ply file name&gt;  
The output file name. It should be in .ply format.

\[\-\-iters &lt;maximum number of iterations&gt;\]  
The maximum number of iterations. The default value of this parameter is 30.

\[\-\-pointWeight &lt;interpolation weight&gt;\]  
The pointWeight parameter of screened Poisson surface reconstruction. The default value for this parameter is 10.

\[\-\-depth &lt;reconstruction depth&gt;\]  
The depth parameter of screened Poisson surface reconstruction. It is the maximum possible depth of the octree. The default value of this parameter is 10.

\[\-\-neighbors &lt;number of neighbors&gt;\]  
The number of the closest sample points to search from every reconstructed triangle face. The suggested range is between 10 and 20. The default value of this parameter is 10.

#### Parameters:
The models in the data folder all could run with default parameters.
