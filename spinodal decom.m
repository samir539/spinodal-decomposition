function driveCH ()
%
% Set some parameters
noiseAmp = 0.0005;
cellSize = 256;
%
% contourPlots is a flag that determines whether contour plots of the
% final compositions are produced.
% 0 ==> Do not produce contour plots
% 1 ==> Produce contour plots
contourPlots = 1;
%
% Define the range of mesh spacings
dxMin = 1.0;
dxMax = 1.0;
dxStep = 0.2;
nSpacings = 1 + round((dxMax-dxMin)/dxStep);
%
% Define the range of time steps
dtMin = 1.0;
dtMax = 1.0;
dtStep = 0.2;
nSteps = 1 + round((dtMax-dtMin)/dtStep);
%
% Define the range of compositions to be considered
cMin = 0.3;
cMax = 0.7;
cStep = 0.05;
nCompositions = 1 + round((cMax-cMin)/cStep);
%
% Define the range of simulation times to be computed
tMin = 100;
tMax = 1000;
tStep = 10;
nTimes = 1 + round((tMax-tMin)/tStep);
%
% Allocate storage
nHeader = 3;
nSets = nSpacings*nSteps*nCompositions;
frac = zeros(nHeader+nTimes,1+nSets);
%
% Allocate times
runTime = tMin;
for i = 1:nTimes
    frac(nHeader+i,1) = runTime;
    runTime = runTime + tStep;
end
%
% Loop over variables and perform simulations
n = 2;
for meshSpacing = dxMin:dxStep:dxMax
    frac(1,n) = meshSpacing;
    for timeStep = dtMin:dtStep:dtMax
        frac(2,n) = timeStep;
        for initialComp = cMin:cStep:cMax
            frac(3,n) = initialComp;
            fprintf ('dx = %6.4f  dt = %6.4f  C0 = %6.4f > ', meshSpacing, timeStep, initialComp);
            for i = 1:nTimes
                runTime = frac(nHeader+i,1);
                finalComp = cahnHilliard(initialComp, noiseAmp, cellSize, meshSpacing, timeStep, runTime);
                frac(nHeader+i,n) = countSites (finalComp);
                fprintf ('.');
            end
            fprintf ('\n');
            if contourPlots == 1
                figure;
                contourf(finalComp);
                title(['dx = ', num2str(meshSpacing), ' dt = ', num2str(timeStep), ' C0 = ', num2str(initialComp)]);
            end
            n = n + 1;
        end
    end
end
%
% Write out results for plotting
writeData(nSets, nTimes, nHeader, frac);
end

function frac = countSites (comp)
[n, m] = size(comp);
nPhase1 = 0;
nPhase2 = 0;
for i = 1:n
    for j = 1:m
        if comp(i,j) < 0.2
            nPhase1 = nPhase1 + 1;
        elseif comp(i,j) > 0.8
            nPhase2 = nPhase2 + 1;
        end
    end
end
frac = (nPhase1+nPhase2)/(n*m);
end

function writeData (nSets, nTimes, nHeader, frac)
%
% Open file to save the statistics
file = fopen ('results.csv', 'w');
%
% Write the header
fprintf(file, 'Mesh spacing');
for i = 1:nSets
    fprintf(file, ', dx = %8.4f', frac(1,1+i));
end
fprintf(file, '\n');
%
fprintf(file, 'Time step');
for i = 1:nSets
    fprintf(file, ', dt = %8.4f', frac(2,1+i));
end
fprintf(file, '\n');
%
fprintf(file, 'Composition');
for i = 1:nSets
    fprintf(file, ', C0 = %8.4f', frac(3,1+i));
end
fprintf(file, '\n');
%
% Write the results
for j = 1:nTimes
    fprintf(file, '%8.4f', frac(nHeader+j,1));
    for i = 1:nSets
        fprintf(file, ' ,%8.4f', frac(nHeader+j,i+1));
    end
    fprintf(file, '\n');
end
%
% Close the file
fclose(file);
end

function finalComp = cahnHilliard( initialComp, noiseAmp, cellSize, meshSpacing, timeStep, runTime )
%CAHNHILLIARD Performs a Cahn-Hilliard simulation of spinidal decomposition
%   The Cahn-Hilliard equations allow spinodal decomposition to be modelled
%   by means of a partial differential equation. The algorithm used here is
%   the one implemented in the prgram CH-muse, written by Gururajan
%   Mogadalai. See:
%   http://sites.google.com/site/gurusofficialhomepage/research/downloads/write-ups-and-codes/the-art-of-phase-field-modelling/ch-muse
%
%   Input arguments are:
%   initialComp       A number between 0 and 1 giving the average composition
%                     of the sample
%   noiseAmp          The amplitude of the random noise added to the initial
%                     composition
%   cellSize          Length of the sides of the simulation box
%   meshSpacing       The spacing of the mesh points used to represent the
%                     composition
%   timeStep          The time step for the integrator
%   runTime           The duration of the simulation
%
%   Output variables:
%   finalComp         The final composition distribution
%
%   Units used in this program:
%   Length and time are in units defined by the coefficients in the
%   Cahn-Hilliard equation. Thus the equations here are dimensionless.

% Get the mesh size for the noise added to the initial
% concentration. Note that we insist on the number of mesh points being an
% even number to avoid problems with fast fourier transforms.
% noiseMeshSpacing is the spacing of the noise mesh points. This is turn is
% a measure of the highest spatial frequency the noise can have.
noiseMeshSpacing = 5.0;
noiseNMesh = max(1,round(cellSize/noiseMeshSpacing));
noiseNMesh = noiseNMesh + mod(noiseNMesh, 2);

% Get the number of mesh points for the solver mesh. We require it
% to be a multiple of the number of points on the noise mesh so
% that we can map the course noise data onto the finer solver mesh.
ratio = max(1,round(noiseMeshSpacing/meshSpacing));
nMesh = noiseNMesh*ratio;

% Build the initial concentration
comp = initialComp*ones(nMesh);
for i = 0:noiseNMesh-1
    ic1 = 1+i*ratio;
    for j = 0:noiseNMesh-1
        jc1 = 1+j*ratio;
        comp(ic1:ic1+ratio-1, jc1:jc1+ratio-1) = (initialComp+normrnd(0.0, noiseAmp))*ones(ratio);
    end
end

% Build the list of g vectors
g1 = fftFreq(nMesh, cellSize);
g1 = g1.*g1;
g2 = zeros(nMesh);
for i=1:nMesh
    for j = 1:nMesh
        g2(i,j) = g1(i) + g1(j);
    end
end
g4 = g2.*g2;

% Propagate the equation of motion
nSteps = round(runTime/timeStep) + 1;
ftComp = fft2(comp);
for i=1:nSteps
    h = comp.*(1.0-comp).*(1.0-2.0*comp);
    ftH = fft2(h);
    ftComp = (ftComp-timeStep*g2.*ftH)./(1.0+timeStep*g4);
    comp = real(ifft2(ftComp));
end
finalComp = comp;
end

function omega = fftFreq( N, extent )
%FFTFREQ returns a vector of the angular frequencies for an FFT
%
%   Input arguments are:
%   N       The number of sample points
%   extent  The full duration of the original sampled signal
%
%   Output arguments are:
%   omega   The vector of angular frequencies
%
omega0 = 2.0*pi/extent;
omega1 = N*omega0;
omega = [0:N-1]*omega0;
for i=N/2+1:N
    omega(i) = omega(i) - omega1;
end
end
%adapted from Prof. Horsfield's course on Materials Modelling