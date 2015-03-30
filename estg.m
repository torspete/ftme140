function G=estg(ptype,e)
%ESTG: Establish the matrix G from Euler angles or Euler parameters
%Inputs: ptype - Parameter type
%        e     - The 3-element Euler angle vector OR the 4-parameter
%                Euler parameter vector
%Output: G     - The G matrix
%Call:   G=estg(ptype,e)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2015-02-21

switch lower(ptype)
  case 'eulera'
    fi=e(1);th=e(2);%ps=ea(3);
    cfi=cos(fi);sfi=sin(fi);cth=cos(th);sth=sin(th);
    G=[0 cfi  sth*sfi;
       0 sfi -sth*cfi;
       1  0     cth];
  case 'eulerp'
    G=2*este(e);
  otherwise
    error
end    
 
   