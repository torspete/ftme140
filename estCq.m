function [Cq,Cqd,Qt]=estCq(z,const,body)
%% ESTCQ: Establish the Jacobian of the constraints
%Inputs: z      - The extended state vector
%        const  - The structured array const
%        body   - The structured array body
%Output: Cq     - The partial derivative of C wrt q
%        Cqd    - The partial derivative of C wrt qdot (equals the partial
%                 derivative of the "hidden" constraint wrt q)
%        Qt     - The partial derivative of C wrt time (with negative sign)
%Call:   [Cq,Cqd,Qt]=estCq(z,const,body)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-03-27
%Modified: 2009-05-19 (Corrected Cq for revolute joint)
%Modified: 2009-07-01 (Corrected Cq and Cqd for translational, sliding and
%                      universal joints)

%%                                                                 Initiate
nb=body.nb;
nceq=const.nceq;
nj=const.nj;
nq=7*nb;
Cq=zeros(nceq,nq);Cqd=Cq;Qt=sparse(zeros(2*nq+nceq,1));

%%                             Establish the matrices for joint constraints
% ncdone=nb;
ncdone=0;
for J=1:nj
  ctype=const.type{J};
  body1=const.body1{J};
  body2=const.body2{J};
  dofs1=(body1-1)*7+[1:7];
  dofs2=(body2-1)*7+[1:7];
  if strcmpi(ctype(1:3),'FIX')        
      ei =z((body1-1)*7+[4:7]);Ei=este(ei);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);
      rib=const.r1{J};    
      if body2>0
        ej =z((body2-1)*7+[4:7]);ejd=z((nb+body2-1)*7+[4:7]);
      else
        ej=[1 0 0 0]';ejd=[0 0 0 0 ]';
      end
      Ej=este(ej);Ejd=este(ejd);      
      Cq(ncdone+[1:6],dofs1)= ...
          [zeros(3,3) -2*Eid*tilde4(rib);
           zeros(3,3) -2*Eid            ];
      Cqd(ncdone+[1:6],dofs1)= ...
          [eye(3,3) -2*Ei*tilde4(rib) ;
           zeros(3,3) 2*Ei            ];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:6],dofs2)= ...
          [zeros(3,3)  2*Ejd*tilde4(rjb) ;
           zeros(3,3)  2*Ejd             ];          
        Cqd(ncdone+[1:6],dofs2)= ...
          [-eye(3,3)  2*Ej*tilde4(rjb) ;
           zeros(3,3) -2*Ej            ];          
      end
      ncdone=ncdone+6;
  elseif strcmpi(ctype(1:3),'REV')        
      ei =z((body1-1)*7+[4:7]);Ei=este(ei);Eib=este(ei,1);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);Eibd=este(eid,1);
      eic=const.e1{J};Aic=esta('eulerp',eic);Aic=Aic(:,1:2);
      rib=const.r1{J};
      if body2>0
        ej =z((body2-1)*7+[4:7]);
        ejd=z((nb+body2-1)*7+[4:7]);
      else
        ej=[1 0 0 0]';
        ejd=[0 0 0 0 ]';
      end
      Ej=este(ej);Ejd=este(ejd);
      
      Cq(ncdone+[1:5],dofs1)= ...
          [zeros(3,3)  -2*Eid*tilde4(rib)   ;
           zeros(2,3) -2*Aic'*(Eibd+2*Eib*tilde4(Ej*ejd,1))];% Corrected
      Cqd(ncdone+[1:5],dofs1)= ...
          [eye(3,3) -2*Ei*tilde4(rib)  ;
          zeros(2,3) 2*Aic'*Eib];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:5],dofs2)= ...
          [zeros(3,3)  2*Ejd*tilde4(rjb)   ;
           zeros(2,3)  2*Aic'*Ei*Eib'*Ejd];          
        Cqd(ncdone+[1:5],dofs2)= ...
          [-eye(3,3)  2*Ej*tilde4(rjb)   ;
          zeros(2,3) -2*Aic'*Ei*Eib'*Ej];          
      end
      ncdone=ncdone+5;
      
  elseif strcmpi(ctype(1:3),'TRA')
      Rid=z((nb+body1-1)*7+[1:3]);
      ei =z((body1-1)*7+[4:7]);Ai=esta('eulerp',ei);Ei=este(ei);Eib=este(ei,1);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);
      eic=const.e1{J};Aic=esta('eulerp',eic);Aic=Aic(:,1:2);
      ej =z((body2-1)*7+[4:7]);Aj=esta('eulerp',ej);Ejb=este(ej,1);
      ejc=const.e2{J};Ajc=esta('eulerp',ejc);Ajc=Ajc(:,1:2);
      rib=const.r1{J};rjb=const.r2{J};
      if body2>0
        Rjd=z((nb+body2-1)*7+[1:3]);
        ej =z((body2-1)*7+[4:7]);
        ejd=z((nb+body2-1)*7+[4:7]);
      else
        Rjd=zeros(3,1);
        ej=[1 0 0 0]';
        ejd=zeros(4,1);
      end
      Ej=este(ej);Ejd=este(ejd);
      
      Cq(ncdone+[1:5],dofs1)= ...
          [zeros(2,3)  Aic'*(Eib*tilde4(Rid,1)+Omega34(eid,rib)) ;
           zeros(3,3)            -2*Eid               ];
      Cqd(ncdone+[1:5],dofs1)= ...
          [Aic'*[Eib*Ei' -2*Ai'*Ei*tilde4(rib)] ;
            zeros(3,3)           2*Ei          ];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:5],dofs2)= ...
          [zeros(2,3)  Ajc'*(Ejb*tilde4(Rjd,1)+Omega34(ejd,rjb)) ;
           zeros(3,3)  2*Ejd                          ];          
        Cqd(ncdone+[1:5],dofs2)= ...
          [-Ajc'*[Ejb*Ej' -2*Aj'*Ej*tilde4(rjb)];
             zeros(3,3)         -2*Ej          ];          
      end
      ncdone=ncdone+5;
  elseif strcmpi(ctype(1:3),'SPH')        
      ei =z((body1-1)*7+[4:7]);Ei=este(ei);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);
      rib=const.r1{J};    
      if body2>0
        ej =z((body2-1)*7+[4:7]);ejd=z((nb+body2-1)*7+[4:7]);
      else
        ej=[1 0 0 0]';ejd=[0 0 0 0 ]';
      end
      Ej=este(ej);Ejd=este(ejd);      
      Cq(ncdone+[1:3],dofs1)= ...
          [zeros(3,3) -2*Eid*tilde4(rib)];
      Cqd(ncdone+[1:3],dofs1)= ...
          [eye(3,3) -2*Ei*tilde4(rib)];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:3],dofs2)= ...
          [zeros(3,3)  2*Ejd*tilde4(rjb)];
        Cqd(ncdone+[1:3],dofs2)= ...
          [-eye(3,3)  2*Ej*tilde4(rjb)];
      end
      ncdone=ncdone+3;
      
  elseif strcmpi(ctype(1:3),'SLI')
      Rid=z((nb+body1-1)*7+[1:3]);
      ei =z((body1-1)*7+[4:7]);Ei=este(ei);Eib=este(ei,1);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);Eibd=este(eid,1);
      eic=const.e1{J};Aic=esta('eulerp',eic);Aic12=Aic(:,1:2);Aic3=Aic(:,3);
      ejc=const.e2{J};Ajc=esta('eulerp',ejc);Ajc3=Ajc(:,3);
      rib=const.r1{J};rjb=const.r2{J};
      if body2>0
        Rjd=z((nb+body2-1)*7+[1:3]);
        ej =z((body2-1)*7+[4:7]);
        ejd=z((nb+body2-1)*7+[4:7]);
      else
        Rjd=zeros(3,1);
        ej=[1 0 0 0]';
        ejd=zeros(4,1);
      end
      Ej=este(ej);Ejd=este(ejd);
      
      Cq(ncdone+[1:3],dofs1)= ...
          [zeros(1,3)  Aic3'*(Eib*tilde4(Rid,1)+Omega34(eid,rib)) ;
           zeros(2,3) -2*Aic12'*(Eibd+Eib*tilde4(Ej*ejd,1))];
      Cqd(ncdone+[1:3],dofs1)= ...
          [Aic3'*[Ai' -2*Ai'*Ei*tilde4(rib)] ;
                zeros(2,3) 2*Aic12'*Eib];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:3],dofs2)= ...
          [zeros(1,3)  Ajc3'*(Ejb*tilde4(Rjd,1)+Omega34(ejd,rjb)) ;
                     zeros(2,3)  2*Aic12'*Ai'*Ejd];          
        Cqd(ncdone+[1:3],dofs2)= ...
          [-Ajc3'*[Aj' -2*Aj'*Ej*tilde4(rjb)];
             zeros(2,3) -2*Aic12'*Ai'*Ej];          
      end
      ncdone=ncdone+3;
  elseif strcmpi(ctype(1:3),'UNI')        
      ei =z((body1-1)*7+[4:7]);Ei=este(ei);Eib=este(ei,1);
      eid=z((nb+body1-1)*7+[4:7]);Eid=este(eid);Eibd=este(eid,1);
      eic=const.e1{J};Aic=esta('eulerp',eic);Aic=Aic(:,3);
      ejc=const.e2{J};Ajc=esta('eulerp',ejc);Ajc=Ajc(:,3);
      rib=const.r1{J};    
      if body2>0
        ej =z((body2-1)*7+[4:7]);ejd=z((nb+body2-1)*7+[4:7]);
      else
        ej=[1 0 0 0]';ejd=[0 0 0 0 ]';
      end
      Ej=este(ej);Ejb=este(ej,1);Ejd=este(ejd);Ejbd=este(ejd,1);      
      Cq(ncdone+[1:4],dofs1)= ...
          [zeros(3,3) -2*Eid*tilde4(rib);
           zeros(1,3) -2*Aic'*Eibd     ];
      Cqd(ncdone+[1:4],dofs1)= ...
          [eye(3,3) -2*Ei*tilde4(rib) ;
           zeros(1,3) 2*Aic'*Eib     ];
      if body2>0
        rjb=const.r2{J};
        Cq(ncdone+[1:4],dofs2)= ...
          [zeros(3,3)  2*Ejd*tilde4(rjb) ;
           zeros(1,3)  2*Ajc'*Ejbd      ];          
        Cqd(ncdone+[1:4],dofs2)= ...
          [-eye(3,3)  2*Ej*tilde4(rjb) ;
           zeros(1,3) -2*Ajc'*Ejb     ];          
      end
      ncdone=ncdone+4;
            
  end    
end

%%                                              Euler parameter constraints
for J=1:nb
  dofs=(J-1)*7+[4:7];
  Cqd(J+ncdone,dofs)=2*z((J-1)*7+[4:7])';
  Cq(J+ncdone,dofs)=2*z((nb+J-1)*7+[4:7])';
end

%%
Cq=sparse(Cq);
Cqd=sparse(Cqd);

function O34=Omega34(e,r)
%% Omega34: The O34 operator for 3D and 4D vectors to 3x4 matrix
%Inputs: e    - A 4-element vector
%        r    - A 3-element vector
%Output: O34  - O34 3x4 matrix
%Call:   Omega34=O34(e,r)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-04-06

e0=e(1);e1=e(2);e2=e(3);e3=e(4);
r1=r(1);r2=r(2);r3=r(3);

O34=-2*[ (+e3*r2-e2*r3) (-e2*r2-e3*r3) (+e1*r2+e0*r3) (-e0*r2+e1*r3);
         (-e3*r1+e1*r3) (+e2*r1-e0*r3) (-e1*r1-e3*r3) (+e0*r1+e2*r3);
         (+e2*r1-e1*r2) (+e3*r1+e0*r2) (-e0*r1+e3*r2) (-e1*r1-e2*r2)];
    
