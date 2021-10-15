function [Value,varargout] = parsevarargin(IN,theoption,thedefault)
% small handmade function for sorting varargin input in one line
% O.Codol 4th May 2019
% codol.olivier@gmail.com
%-------------

k = find(strcmpi(IN,theoption));
if any(k)
  Value = IN{k(1)+1};
  IN( [ k(1) k(1)+1 ] ) = [];
else
    Value = thedefault;
end

if strcmpi(Value,'init'); figure; Value=subplot(1,1,1); end

varargout = IN;