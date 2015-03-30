function [e,v,theta]=a2eulerp(A)
%A2EULERP: Obtain the Euler parameters related to a given rotation matrix A
%Inputs: A   - The rotation matrix A
%Output: e   - The 4-element Euler parameter vector
%Call:   [e,v,theta]=a2eulerp(A)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-01-30

% theta=acos((trace(A)-1)/2);
% 
% % --------------------------------------------------------------------------
% %                                             If rotation is extremely small
% %                                             v does not matter much, set it
% %                                             to [1 0 0]
% %                                             ------------------------------
% if abs(theta)>eps;
%   v=[A(3,2)-A(2,3) A(1,3)-A(3,1) A(2,1)-A(1,2)]/2/sin(theta);
% else
%   v=[1 0 0];
% end
% 
% % --------------------------------------------------------------------------
% %                                                       The Euler parameters
% %                                                       --------------------
% e=[cos(theta/2);sin(theta/2)*v(:)];
    

% --------------------------------------------------------------------------
%                                                       The Euler parameters
%                                                       --------------------
trA=trace(A);
e0=sqrt((trA+1)/4);
e1=sqrt((1+2*A(1,1)-trA)/4);
e2=sqrt((1+2*A(2,2)-trA)/4);
e3=sqrt((1+2*A(3,3)-trA)/4);
if e0>eps
  e1=(A(3,2)-A(2,3))/4/e0;
  e2=(A(1,3)-A(3,1))/4/e0;
  e3=(A(2,1)-A(1,2))/4/e0;
else
  if e1>eps
    e2=(A(2,1)+A(1,2))/4/e1;
    e3=(A(1,3)+A(3,1))/4/e1;
  elseif e2>eps
    e1=(A(2,1)+A(1,2))/4/e2;
    e3=(A(3,2)+A(2,3))/4/e2;
  else
    e1=(A(1,3)+A(3,1))/4/e3;
    e2=(A(3,2)+A(2,3))/4/e3;
  end    
end    
e=[e0 e1 e2 e3]';  