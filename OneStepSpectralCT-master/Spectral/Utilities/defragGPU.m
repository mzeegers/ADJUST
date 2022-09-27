function res = defragGPU(listofvars)
%DEFRAGGPU Defragments memory on the GPU
%   Offloads all variables from the GPU, resets the GPU, then loads the
%   data back
res = 0;

% Get the current GPU device
g = gpuDevice();

% Get the list of variables in memory
listonCPU = {};
listofnames = {};

counter = 0;
% Scroll through the variables
for c = 1:numel(listofvars)
    currentvar= evalin('caller', listofvars{c});
    % If it is a gpuArray
    if isa(currentvar, 'gpuArray')
        counter = counter + 1;
        % Offload it from the GPU
        listonCPU{counter} = gather(currentvar);
        listofnames{counter} = listofvars{c};
    end
end

% Reset the current GPU device
reset(g);

% Load the data back to the GPU device
for c = 1:numel(listonCPU)
    assignin('caller', listofnames{c}, gpuArray(listonCPU{c}));
end

end

