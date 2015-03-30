function [stri,inte,real,comp,name]=decodestr(line,delim)
% DECODESTR convert a line of combined text and numbers into components.
%
% [stri,inte,real,comp,name]=decodestr(line,delim)
% Input:    line - string to be deassembled
%           delim - (optional) delimiter
% Output:   stri - vector of character strings in order of appearance
%                  (maximum 16 character strings)
%           inte - vector of integers in order of appearance
%           real - vector of reals in order of appearance
%           comp - vector of complex numbers in order of appearance
%           name - The name(s) given within ':s on line, e.g. 'Thomas'

% Copyright: Thomas Abrahamsson, Chalmers University of Technology, Sweden
% Written:  March 25, 2009

if nargin>1
  if length(delim)>1,error('Delimiter needs to be single character');end
  if ~ischar(delim)>1,error('Delimiter needs to be single character');end  
  ind=findstr(line,delim);line(ind)=' ';
end

% -----------------------------------------------------------------------------
%                                                1. Add blanks to line
%                                                ------------------------------
line=[' ' line ' '];

% -----------------------------------------------------------------------------
%                                                2. Identify all blanks
%                                                ------------------------------
bl=find(line==' ');

% -----------------------------------------------------------------------------
%                                                3. Remove excessive blanks
%                                                ------------------------------
addr=find(line==' ');line(addr(find(diff([-1 addr])==1)))=[];

% -----------------------------------------------------------------------------
%                                                4. Identify all blanks again
%                                                ------------------------------
bl=find(line==' ');

% -----------------------------------------------------------------------------
%                                                5. Identify components
%                                                ------------------------------
stri=[];inte=[];real=[];comp=[];
for I=1:length(bl)-1,
  test_string=line(bl(I)+1:bl(I+1)-1);
  num=str2double(test_string);if isnan(num),num=[];end
  lstr=length(test_string);
  if isempty(num)
    tmp='                ';
    tmp(1:min([lstr 16]))=test_string(1:min([lstr 16]));
    stri=[stri;tmp];
  elseif (imag(num))~=0,
    comp=[comp;num];
  else
    if any(test_string=='.') | any(lower(test_string)=='e'),
      real=[real;num];
    else
      inte=[inte;num];
    end
  end
end

% -----------------------------------------------------------------------------
%                                            5. Check if any strings were names
%                                            ----------------------------------
name=[];
for I=size(stri,1):-1:1
    str=deblank(stri(I,:));
    if strcmp(str(1),'''') & strcmp(str(end),'''')
       name=char(name,str(2:end-1));
       stri(I,:)=[];
   end
end    
if ~isempty(name),name(1,:)=[];end