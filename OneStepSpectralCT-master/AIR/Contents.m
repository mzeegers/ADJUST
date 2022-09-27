% AIR Tools.
% Version 1.1  11-August-13
%
% Iterative ART Methods.
%   kaczmarz        - Kaczmarz's method (often referred to as ART).
%   randkaczmarz    - Randomized Kaczmarz method.
%   symkaczmarz     - Symmetric Kaczmarz method.
%
% Iterative SIRT Methods.
%   cav             - Component Averaging (CAV) method.
%   cimmino         - Cimmino's method.
%   drop            - Diagonally Relaxed Orthogonal Projections (DROP) method.
%   landweber       - The classical Landweber method.
%   sart            - Simultaneous Algebraic Reconstruction Technique (SART) method.
%
% Training Routines.
%   trainDPME       - Training method for the stopping rules DP and ME.
%   trainLambdaART  - Training to determine optimal lambda for ART methods.
%   trainLambdaSIRT - Training to determine optimal lambda for SIRT method.
%
% Test Problems.
%   binarytomo      - 2D tomography test problem with a binary image.
%   fanbeamtomo     - 2D tomography test problem using fan beams.
%   odftomo         - 2D tomography test problem with a smooth odf image.
%   paralleltomo    - 2D tomography test problem using parallel beams. 
%   seismictomo     - 2D seismic tarvel-time tomography test problem.
%
% Demonstration Scripts.
%   ARTdemo         - Demonstrates the use of, and the results from, the ART methods.
%   nonnegdemo      - Demonstrates the use of nonnegativity constraints.
%   SIRTdemo        - Demonstrates the use of, and the results from, the SIRT methods.
%   trainingdemo    - Demonstrates the use of the training methods.
%
% Auxiliary Routines.
%   calczeta        - Calculates a specific root of a certain polynomial.
%   rzr             - Remove zero rows of A and the corresponding elements of b.