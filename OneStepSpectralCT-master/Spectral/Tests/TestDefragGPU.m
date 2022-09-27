function [sizeA] = TestDefragGPU()
A = gpuArray(magic(100));
defragGPU(who());
sizeA = size(A);
end

