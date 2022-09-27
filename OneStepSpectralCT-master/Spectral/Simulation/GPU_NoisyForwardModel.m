function [y, k_d] = GPU_NoisyForwardModel(object, A, M, S, noisy, normalize)

% Convert all inputs to GPU arrays
object = gpuArray(object);
A = gpuArray(A);
M = gpuArray(M);
S = gpuArray(S);

% Compute forward projections
projs = A * object;

% Compute photon counts as seen by the detector
attenuationFactors = exp(- M * projs.');
counts = gather(S * attenuationFactors);
counts_without_attenuation = gather(S * ones(size(attenuationFactors), 'gpuArray'));


if (noisy)
	% Add noise to the photon counts
	noisy_counts = poissrnd(counts);
else
	noisy_counts = counts;
end

if (normalize)
	% Divide by the counts in the absence of object
	y = noisy_counts ./ counts_without_attenuation;
else
	y = noisy_counts;
end

y_without_attenuation = poissrnd(counts_without_attenuation) ./ counts_without_attenuation;
k_d = var(y_without_attenuation(:)) / mean(y_without_attenuation(:));

end
