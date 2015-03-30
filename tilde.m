function vt=tilde(v)
%VT: The tilde operator for 3D vector to 3x3 matrix
%Inputs: v   - A 3-element vector
%Output: vt  - The tilde 3x3 matrix
%Call:   vt=tilde(v)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-01-30

try
  vt=sparse([3;1;2],[2;3;1],v,3,3);vt=full(vt-vt');
catch
    error(lasterr)
end 
