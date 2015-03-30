function vt4=tilde4(v,opt)
%VT4: The tilde4 operator for 3D vector to 4x4 matrix
%Inputs: v   - A 3-element vector
%        opt - Option flag, if opt=1 then tilde-tilde-4 operation is made
%Output: vt4 - The tilde4 4x4 matrix
%Call:   vt4=tilde4(v,opt)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-27

if nargin<2,opt=0;end

if opt==0
  try
    vt4=[  0     v(:)'  ;
       -v(:) tilde(v) ];
  catch
     error(lasterr)
  end
elseif opt==1
  try
    vt4=[  0     v(:)'  ;
       -v(:) -tilde(v) ];
  catch
     error(lasterr)
  end   
end