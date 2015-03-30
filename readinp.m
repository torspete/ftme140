function [Body,Const,Spring,Case]=readinp(infile)
%READINP: Read the input file
% See helpmubody (help helpmubody) for description of objects
% Body, Const, Spring and Case

% Copyleft: Thomas Abrahamsson, Chalmers University of Technology, Sweden
% Written: 2009-03-24 /TA
% Modified: 2009-03-26 /TA (added option to include general spring stiffness)
% Modified: 2009-04-06 /TA (added option to include initial rates)
% Modified: 2009-04-07 /TA (corrected errors)

% -------------------------------------------------------------------------
%                                                                  Initiate
%                                                                  --------
warnstate=warning;warning('off');
Disctimes=[];
Body.nb=0;Const.nj=0;nj=0;Spring.ns=0;

% -------------------------------------------------------------------------
%                                                     Open file for reading
%                                                     and read all inputs
%                                                     ---------------------
fid=fopen(infile,'r'); if fid<0,error('File does not exist');end
firstline=1;readnewline=1;
while 1
  if readnewline,line = myfgetl(fid);end,readnewline=1;
  [Stri,Inte,Real,Comp,Name]=decodestr(line);if size(Stri,1)<5,Stri(5,1)=' ';end
  if firstline
     firstline=0;
  elseif strcmp(line(1),'*'),% Do nothing, it's just a comment
  elseif isempty(deblank(line)),% Do nothing, it's just an empty line    

% -------------------------------------------------------------------------
%                                                     BODY MASS AND INERTIA
%                                                     ---------------------   
  elseif strcmpi(Stri(1,1),'B') & strcmpi(Stri(2,1),'M')
     while 1,
       try, nb=length(Body.m);catch,nb=0;end
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else,break,end,end
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if isempty(Stri) & length(Real)==7
         if isempty(Name),Body.name{nb+1}=['Body_' int2str(nb+1)];else,Body.name{nb+1}=Name(1,:);end
         Body.m{nb+1}=Real(1);
         Body.I{nb+1}=[Real(2) Real(3) Real(4);
                       Real(3) Real(5) Real(6);
                       Real(4) Real(6) Real(7)];                   
       else
         break
       end
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else,break,end,end
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if size(Stri,1)==2,Stri(3,1)='R';end
       if size(Stri,1)==3 & length(Real)>5 & length(Real)<8
          Body.R0{nb+1}=[Real(1) Real(2) Real(3)]';
          if upper(Stri(1))=='R' & upper(Stri(2))=='V' & upper(Stri(3))=='R'
              A=esta('rotvect',Real(4:7));
          elseif upper(Stri(1))=='R' & upper(Stri(2))=='V' & upper(Stri(3))=='D'
              A=esta('rotvect',[Real(4:6);pi*Real(7)/180]);
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='P'
              A=esta('eulerp',Real(4:7));
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='A' & upper(Stri(3))=='R'
              A=esta('eulera',Real(4:6));
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='A' & upper(Stri(3))=='D'
              A=esta('eulera',pi*Real(4:6)/180);
          elseif upper(Stri(1))=='R' & upper(Stri(2))=='P'
              A=esta('rodriguez',Real(4:6));
          else
              error('Syntax error on BODY MASS AND INERTIA input card')
          end
          Body.e0{nb+1}=a2eulerp(A);      
       else
         error('Syntax error on BODY MASS AND INERTIA input card')
       end
       Body.nb=Body.nb+1;
     end
     Const.nceq=Body.nb;
         
% -------------------------------------------------------------------------
%                                                         JOINT CONSTRAINTS
%                                                         -----------------   
  elseif strcmpi(Stri(1,1),'J') & strcmpi(Stri(2,1),'C')
    try Const.nceq; catch error('Bodies need to be specified before Joint Constraints');end
    while 1,
      %try nc=length(Const.body1);catch,nc=0;end
      line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else break,end,end
      readnewline=0;
      [Stri,Inte,Real,Comp,Name]=decodestr(line);
      if size(Stri,1)<4,Stri(4,1)='R';end
      if size(Stri,1)==4 & length(Inte)==1 & length(Real)>2
        nj=nj+1;
        Const.nj=Const.nj+1;
        if isempty(Name),Const.name{nj}=['Joint_' int2str(nj+1)];else Const.name{nj}=Name(1,:);end
        Const.body1{nj}=Inte(1);
        Const.r1{nj}=Real(1:3);
        JointType=upper(Stri(1,1:3));Const.type{nj}=JointType;
        if strcmp(JointType,'FIX') | strcmp(JointType,'SPH') 
          Const.e1{nj}=[1 0 0 0];
        else
          if upper(Stri(2))=='R' & upper(Stri(3))=='V' & upper(Stri(4))=='R'
            A=esta('rotvect',Real(4:7));
          elseif upper(Stri(2))=='R' & upper(Stri(3))=='V' & upper(Stri(4))=='D'
            A=esta('rotvect',[Real(4:6);pi*Real(7)/180]);
          elseif upper(Stri(2))=='E' & upper(Stri(3))=='P'
            A=esta('eulerp',Real(4:7));
          elseif upper(Stri(2))=='E' & upper(Stri(3))=='A' & upper(Stri(4))=='R'
            A=esta('eulera',Real(4:6));
          elseif upper(Stri(2))=='E' & upper(Stri(3))=='A' & upper(Stri(4))=='D'
            A=esta('eulera',pi*Real(4:6)/180);
          elseif upper(Stri(2))=='R' & upper(Stri(3))=='P'
            A=esta('rodriguez',Real(4:6));
          else
            error('Syntax error on JOINT CONSTRAINT input card')
          end
          Const.e1{nj}=a2eulerp(A);          
        end  
      else
        break
      end
      line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else,break,end,end
      readnewline=0;
      [Stri,Inte,Real,Comp,Name]=decodestr(line);
      if size(Stri,1)<3,Stri(3,1)='R';end
      if length(Inte)==1
        Const.body2{nj}=Inte(1);
        if strcmp(JointType,'UNI'),%Orientation at 2nd body only needed for universal joint
          if upper(Stri(1))=='R' & upper(Stri(2))=='V' & upper(Stri(3))=='R'
            A=esta('rotvect',Real(1:4));
          elseif upper(Stri(1))=='R' & upper(Stri(2))=='V' & upper(Stri(3))=='D'
            A=esta('rotvect',[Real(1:3);pi*Real(4)/180]);
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='P'
            A=esta('eulerp',Real(1:4));
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='A' & upper(Stri(3))=='R'
            A=esta('eulera',Real(1:3));
          elseif upper(Stri(1))=='E' & upper(Stri(2))=='A' & upper(Stri(3))=='D'
            A=esta('eulera',pi*Real(1:3)/180);
          elseif upper(Stri(1))=='R' & upper(Stri(2))=='P'
            A=esta('rodriguez',Real(1:3));
          else
            error('Syntax error on JOINT CONSTRAINT input card')
          end
          Const.e2{nj}=a2eulerp(A);                     
        end
      else
        error('Syntax error on JOINT CONSTRAINT input card')
      end
      if strcmp(JointType,'FIX')
          Const.ind{nj}=Const.nceq+[1:6];Const.nceq=Const.nceq+6;
      elseif strcmp(JointType,'TRA')
          Const.ind{nj}=Const.nceq+[1:5];Const.nceq=Const.nceq+5;
      elseif strcmp(JointType,'REV')
          Const.ind{nj}=Const.nceq+[1:5];Const.nceq=Const.nceq+5;
      elseif strcmp(JointType,'SPH')
          Const.ind{nj}=Const.nceq+[1:3];Const.nceq=Const.nceq+3;
      elseif strcmp(JointType,'SLI')
          Const.ind{nj}=Const.nceq+[1:4];Const.nceq=Const.nceq+4;
      elseif strcmp(JointType,'UNI')
          Const.ind{nj}=Const.nceq+[1:5];Const.nceq=Const.nceq+5;
      end    
    end

% -------------------------------------------------------------------------
%                                                           SPRING ELEMENTS
%                                                           ---------------   
  elseif strcmpi(Stri(1,1),'S') & strcmpi(Stri(2,1),'E')
     while 1,
       try, ns=length(Spring.k);catch,ns=0;end
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else break,end,end
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if length(Inte)==2 & length(Real)==8
         if isempty(Name),Spring.name{ns+1}=['Spring_' int2str(ns+1)];else Spring.name{ns+1}=Name(1,:);end
         if size(Name,1)~=2,Spring.Ffun{ns+1}=[];else Spring.Ffun{ns+1}=Name(2,:);end
         Spring.k{ns+1}=Real(1);
         Spring.L0{ns+1}=Real(2);
         Spring.body1{ns+1}=Inte(1);
         Spring.r1{ns+1}=Real(3:5);
         Spring.body2{ns+1}=Inte(2);
         Spring.r2{ns+1}=Real(6:8);
       else
         break
       end
       Spring.ns=Spring.ns+1;
     end
             
% -------------------------------------------------------------------------
%                                                     GRAVITY SPECIFICATION
%                                                     ---------------------   
  elseif strcmpi(Stri(1,1),'G') & strcmpi(Stri(2,1),'S')
     while 1,
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else,break,end,end
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if length(Real)==3
         Case.gravv=Real(1:3);
       else
         break
       end       
     end

% -------------------------------------------------------------------------
%                                               INHOMOGENEOUS INITIAL RATES
%                                               ---------------------------   
  elseif strcmpi(Stri(1,1),'I') & strcmpi(Stri(2,1),'I') & strcmpi(Stri(3,1),'R')
     while 1,       
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else break,end,end
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if length(Inte)==1 & length(Real)==7
          J=Inte(1);
          Body.R0d{J}=Real(1:3);
          Body.e0d{J}=Real(4:7);
       else
         break
       end
     end
     
% -------------------------------------------------------------------------
%                                                       EVALUATION AT TIMES
%                                                       -------------------   
  elseif strcmpi(Stri(1,1),'E') & strcmp(Stri(2,1),'A') & strcmpi(Stri(3,1),'T')
     while 1,
       line=myfgetl(fid);while 1,if line(1)=='*',line=myfgetl(fid);else break,end,end
       line=strrep(lower(line),'step','STE'); 
       readnewline=0;
       [Stri,Inte,Real,Comp,Name]=decodestr(line);
       if isempty(Stri) | all(Stri(:,1)=='S')
          eval('Disctimes=[Disctimes stepp(line)];','Disctimes=stepp(line);')
       else
          Case.evalt=Disctimes;
          Stri(5,1)=' ';
          break
       end
     end
    
% -------------------------------------------------------------------------
%                                                               END OF DATA
%                                                               -----------
  elseif strcmpi(Stri(1,1),'E') & strcmpi(Stri(2,1),'O') & strcmpi(Stri(3,1),'D')
    break          
  else
     error(['Does not recognize: ' line])          
  end
end

% -------------------------------------------------------------------------
%                                                                Close file
%                                                                ----------
fclose(fid);

warning(warnstate);


function line=myfgetl(fid)
line=fgetl(fid);
line(find(line==9))=' ';% Replace all Tabs with Blanks