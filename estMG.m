function [M,G]=estMG(z,body)
%ESTMG: Establish the mass and gyroscopic matrices for bodies
%Inputs: z     - The extended state vector
%        body  - The structured array body
%Output: M     - The mass matrix
%        G     - The gyroscopic matrix
%Call:   [M,G]=estMG(z,body)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-26

% -------------------------------------------------------------------------
%                                                                  Initiate
%                                                                  --------
nb=body.nb;
nceq=length(z)-14*nb;

% -------------------------------------------------------------------------
%                                                     Establish mass matrix
%                                                     ---------------------
for J=1:nb
    dofsR=(J-1)*7+[1:3];
    dofse=(J-1)*7+[4:7];
    Ebar=este(z(dofse),1);
    M(dofsR,dofsR)=body.m{J}*eye(3,3);
    M(dofse,dofse)=4*Ebar.'*body.I{J}*Ebar;
end    

% -------------------------------------------------------------------------
%                                               Establish gyroscopic matrix
%                                               ---------------------------
for J=1:nb
    dofse=(J-1)*7+[4:7];
    dofsed=(nb+J-1)*7+nceq+[4:7];
    Ebar=este(z(dofse),1);    
    Ebard=este(z(dofsed),1);    
    G(dofse,dofse)=8*Ebard.'*body.I{J}*Ebar;
end

M=sparse(M);
G=sparse(G);
