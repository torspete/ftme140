function G=estgba(ea)
%ESTGBA: Establish the matrix G from Bryant angles
%Inputs: ea    - The 3-element Bryant angle vector
%Output: G     - The G matrix
%Call:   G=estgba(ea)

%Written: 2015-02-21, Thomas Abrahamsson, Chalmers University of Technology

fi=ea(1);th=ea(2);
cf=cos(fi);sf=sin(fi);st=sin(th);ct=cos(th);
G=[1  0    st  ;
   0  cf -sf*ct;
   0  sf  cf*ct];
 
   