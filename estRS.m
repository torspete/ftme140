function [R,S]=estRS(M,K,G,Cq,Cqd)
%ESTRS: Establish the state-space matrices R and S
%Inputs: M      - The mass matrix
%        K      - The stiffness matrix
%        G      - The gyroscopic matrix
%        Cq     - The Jacobian of the constraints
%        Cqd    - The partial derivative of C wrt qdot
%Output: R,S    - The state-space matrices
%Call:   [R,S]=estRS(M,K,G,Cq,Cqd)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-27

% -------------------------------------------------------------------------
%                                                                  Initiate
%                                                                  --------
nq=size(M,1);
nceq=size(Cq,1);
R=zeros(2*nq+nceq,2*nq+nceq);S=R;

% -------------------------------------------------------------------------
%                                                    Establish the matrices
%                                                    ----------------------
dofr=1:nq;
dofc=1:nq;R(dofr,dofc)=eye(nq,nq);
dofc=nq+[1:nq];S(dofr,dofc)=-eye(nq,nq);

dofr=nq+[1:nq];
dofc=nq+[1:nq];R(dofr,dofc)=M;
dofc=2*nq+[1:nceq];R(dofr,dofc)=Cqd';
dofc=[1:nq];S(dofr,dofc)=K;
dofc=nq+[1:nq];S(dofr,dofc)=G;

dofr=2*nq+[1:nceq];
dofc=1:nq;R(dofr,dofc)=Cq;
dofc=nq+[1:nq];R(dofr,dofc)=Cqd;

