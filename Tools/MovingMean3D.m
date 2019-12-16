function out = MovingMean3D(input , Kernel)
  if nargin < 2
    Kernel = ones(5,5,2) /50;
  end
  out = convn(input, Kernel, 'same');
end
