function E=este(e,opt)
%ESTE: Establish the matrix E from Euler parameters
%Inputs: e   - The 4-element Euler parameter vector
%        opt - If given as opt=1, The Ebar matrix is established
%Output: E   - The E matrix (default)
%        Eb  - The Ebar matrix (if opt=1) 
%Call:   E=este(e[,opt])

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-01-30

if nargin<2,opt=0;end

e0=e(1);e1=e(2);e2=e(3);e3=e(4);
if opt==0
  E=[ -e1 e0 -e3 e2;
      -e2 e3 e0 -e1;
      -e3 -e2 e1 e0];
elseif opt==1
  E=[ -e1 e0 e3 -e2;
      -e2 -e3 e0 e1;
      -e3 e2 -e1 e0];
else
    error('Unknown option');
end    
