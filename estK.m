function [K,Qs]=estK(z,spring,body,const)
%ESTK: Establish the stiffness matrix for springs
%Inputs: z      - The extended state vector
%        spring - The structured array spring
%        body   - The structured array body
%Output: K      - The stiffness matrix
%        Qs     - The external force caused by springs connected to ground
%Call:   [K,Qs]=estK(z,spring,body,const)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-26

% -------------------------------------------------------------------------
%                                                                  Initiate
%                                                                  --------
ns=spring.ns;
nq=7*body.nb;
nceq=const.nceq;
K=zeros(nq,nq);
Qs=zeros(2*nq+nceq,1);

% -------------------------------------------------------------------------
%                                                Establish stiffness matrix
%                                                --------------------------
for J=1:ns
    body1=spring.body1{J};
    body2=spring.body2{J};
    rib=spring.r1{J};
    rjb=spring.r2{J};
    Ri=z((body1-1)*7+[1:3]);
    if body2==0,Rj=[0 0 0]';else,Rj=z((body2-1)*7+[1:3]);end
    ei=z((body1-1)*7+[4:7]);
    if body2==0,ej=[1 0 0 0]';else,ej=z((body2-1)*7+[4:7]);end
%   disp('estK'),norm(ei)-1
    ri=Ri+esta('eulerp',ei)*rib;
    rj=Rj+esta('eulerp',ej)*rjb;
    L=norm(ri-rj);    
    kappa=spring.k{J}*(1-spring.L0{J}/L);
    alpha=[         eye(3,3)        ;
           -2*tilde4(rib)'*este(ei)';
                   -eye(3,3)        ;
            2*tilde4(rjb)'*este(ej)'];
    beta=[eye(3,3) -este(ei)*tilde4(rib) -eye(3,3) este(ej)*tilde4(rjb)];
    k=kappa*alpha*beta;
    if body2==0
      dofs=(body1-1)*7+[1:7];
      K(dofs,dofs)=K(dofs,dofs)+k(1:7,1:7);
      Qs(nq+dofs,1)=Qs(nq+dofs,1)-k(1:7,11);
    else    
      dofs=[(body1-1)*7+[1:7] (body2-1)*7+[1:7]];
      K(dofs,dofs)=K(dofs,dofs)+k;
    end  
end

K=sparse(K);
Qs=sparse(Qs);


    
