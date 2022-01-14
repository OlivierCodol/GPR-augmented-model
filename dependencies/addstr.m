function D=addstr(D,A,type,force)
% function D=addstruct(D,A,type,force)
% Adds the fields of structure A to the fields of structure D
% adds the field as a row or column depending on type
% 'last' concatenates a 2D array along columns, or an N-D array along the
%   last dimension.
% 'next' concatenates a 2D array along the 3rd dimension, N-D array along
%   the N+1 dimension.
% type = 'row' / 'column' / 'last / next'
%   row: (DEFAULT) add as rows
%   column: add as columns
%   last: add as the last dimension
%   next: add as the next dimension
% 'force': Forces all fields to be added in same row/column
%          This means, if a field is non-existent or shorter, it will
%          be padded with NaNs
%          If a field is non-existent in the added structure, it will also
%          be padded with NaNs;
% -------------------------------------------------------
if isstruct(A)
    Anames = fieldnames(A);
else
    error('First input class must be a structure.');
end
if isstruct(D)
    Dnames = fieldnames(D);
else
    Dnames = [];
end
if (~isstruct(D) || isempty(Dnames))
    D = struct;
    for k = 1:numel(Anames)
        D.(Anames{k}) = [];
    end
    Dnames = Anames;
end

if nargin <3 || isempty(type)
    dim=1;
    type = [];
elseif ischar(type)
    if strcmp(type,'row')
        dim = 1;
    elseif strcmp(type,'column')
        dim = 2;
    elseif strcmp(type,'last')
        dim = length(size(D.(Dnames{1})));
    elseif strcmp(type,'next')
        dim = length(size(D.(Dnames{1})))+1;
    else
        error('unknown option (row/column/last/next)');
    end
elseif isa(type,'double') && type==round(type)
    dim = type;
else
    error('unknown option (row/column/last/next or a dimension number)');
end

if nargin>3 && strcmp(force,'force'); force = 1; else; force = 0; end

[Dsz,Asz] = deal([]);
for k = 1:numel(Anames)
    if ~isfield(D,Anames{k})
        D.(Anames{k}) = [];
    end
    if force
        Fsz = size(D.(Anames{k}));
        maxsz = nan(2,max(numel(Dsz),numel(Fsz)));
        maxsz(1,1:numel(Fsz)) = Fsz;
        maxsz(2,1:numel(Dsz)) = Dsz;
        Dsz = max(maxsz);
        
        Fsz = size(A.(Anames{k}));
        maxsz = nan(2,max(numel(Asz),numel(Fsz)));
        maxsz(1,1:numel(Fsz)) = Fsz;
        maxsz(2,1:numel(Asz)) = Asz;
        Asz = max(maxsz);
    end
end


for k = 1:numel(Dnames)
    if ~isfield(A,Dnames{k}) && force
        A.(Dnames{k}) = [];
        [A.(Dnames{k}),flag] = forcefcn(A.(Dnames{k}),Asz,[],type);
    else
        flag = 0;
    end
    if flag; warning(['field ' Dnames{k} ' padded']); end
end

        
Anames = fieldnames(A);
for k = 1:numel(Anames)
    if force
        [D.(Anames{k}),flag(1)] = forcefcn(D.(Anames{k}),Dsz,[],type);
        [A.(Anames{k}),flag(2)] = forcefcn(A.(Anames{k}),Dsz,dim,type);
    else
        flag = 0;
    end
    
    D.(Anames{k}) = cat(dim,D.(Anames{k})  ,A.(Anames{k}) );
    if any(flag); warning(['field ' Anames{k} ' padded']); end
end
end

function [Fpad,flag] = forcefcn(F,Ssz,dim,type)
Fsz = size(F);
diffsz = nan(2,max(numel(Fsz),numel(Ssz)));
diffsz(1,1:numel(Fsz)) = Fsz;
diffsz(2,1:numel(Ssz)) = Ssz;
diffsz(isnan(diffsz)) = 1;
if strcmp(type,'next')
    padsz = diff(diffsz);
    padsz(padsz<0) = 0;
    Fpad = padarray(F,padsz,nan,'post');
else
    padsz = diff(diffsz);
    padsz(padsz<0) = 0;
    padsz(dim) = 0;
    Fpad = padarray(F,padsz,nan,'post');
end
padsz = diff(diffsz);
padsz(dim) = 0;
if any(padsz~=0); flag=1; else; flag=0; end
end
