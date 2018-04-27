function D = matStiffTen3D(lambda,mu)
%%MATSTIFFTEN3D produces 3D stiffness tensor for isotropic material
%Input:
%   lambda (scalar): first lame parameter
%   mu (scalar): second lame parametr
%Ouput:
%   D (6x6): stiffness tensor

D11 = lambda + 2*mu;
D = [ D11 lambda lambda 0 0 0;
      lambda D11 lambda 0 0 0;
      lambda lambda D11 0 0 0;
      0 0 0 mu 0 0;
      0 0 0 0 mu 0;
      0 0 0 0 0 mu ];

end