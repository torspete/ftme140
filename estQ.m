function Qx=estQ(z,Body,Case)
%ESTQ: Establish the external load vector
%Inputs: z      - The extended state vector
%        Body   - The structured array body
%        Case   - The structured array case
%Output: Qx     - The external load vector
%Call:   Qx=estK(z,Body,Case)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-30

% -------------------------------------------------------------------------
%                                                                  Initiate
%                                                                  --------
nb=Body.nb;
nz=size(z,1);
Qx=zeros(nz,1);

% -------------------------------------------------------------------------
%                                           Establish the load from gravity
%                                           -------------------------------
g=Case.gravv;
for J=1:nb
  m=Body.m{J};
  Qx((nb+J-1)*7+[1:3])=m*g(:);
end

Qx=sparse(Qx);


    
