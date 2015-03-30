function A=esta(ptype,e)
%ESTA: Establish the rotation matrix A from rotation parameters
%Inputs: ptype - ptype='Eulerp' uses Euler parameters. Default
%                     ='Rodriguez' uses Rodriguez parameters
%                     ='Eulera' uses Euler angles
%                     ='Rotvect' uses rotation vector and rotation angle
%        e     - The rotation parameter vector. Default is 4-element 
%                Euler parameter vector. Other parameter vectors can
%                be specified combined with the optional ptype input
%Output: A     - The rotation matrix A
%Call:   A=esta(ptype,e)

%Copyright: Thomas Abrahamsson, Chalmers, Sweden
%Written: 2009-01-30

if isempty(ptype),ptype='Eulerp';end


try
  if strcmpi(ptype(1:6),'eulerp')
%    if abs(norm(e)-1)>10*eps,norm(e),error('Norm of e is not 1');end
    e=e/norm(e);
    A=este(e)*este(e,1)';
  elseif strcmpi(ptype(1:3),'rod')
      A=eye(3)+2/(1+e'*e)*tilde(e)*(eye(3)+tilde(e));
  elseif strcmpi(ptype(1:6),'eulera')
      fi=e(1);th=e(2);ps=e(3);
      cfi=cos(fi);sfi=sin(fi);cth=cos(th);sth=sin(th);cps=cos(ps);sps=sin(ps);
      D1=[cfi sfi 0;-sfi cfi 0;0 0 1];
      D2=[1 0 0;0 cth sth;0 -sth cth];
      D3=[cps sps 0;-sps cps 0;0 0 1];
      A=(D3*D2*D1)';
  elseif strcmpi(ptype(1:6),'rotvec')
    vt=tilde(e(1:3));   %skew symmetric rotation vector
    vt=vt/norm(vt); % it needs to be a unit vector
    theta=e(4); %rotation angle
    A=eye(3)+vt*sin(theta)+2*vt^2*sin(theta/2)^2;
  else
    error('Does not recognize ptype')
  end
catch
    error(['Could not establish rotation matrix because of error' char(10) lasterr]);
end    
