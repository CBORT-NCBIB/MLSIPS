function intensity = tom2Int(varargin)
%TOM2INT converts channel or Stokes data into intensity data
%
% intensity(tom) computes the intensity from the channel tomogram data
% converts the reconstructed tomogram, either as two channels or two Stokes
% vectors into intensity
%
% Usage:
% ii = tom2Int(S1,S2) where [S1,S2] = recstrTom(filename);
% ii = tom2Int(tom) where tom = recstrTom(filename);
%
% Copyright 2012 Martin Villiger - mvilliger@partners.org

% tomogram with two channels
if nargin==1 && isstruct(varargin{1})
    tom = varargin{1};
    if isfield(tom,'ch1')
        intensity = abs(tom.ch1).^2 + abs(tom.ch2).^2;
    elseif isfield(tom,'ch11')
        intensity = abs(tom.ch11).^2 + abs(tom.ch22).^2 + abs(tom.ch21) + abs(tom.ch21).^2;
%        intensity = abs(tom.ch11.*tom.ch22 - tom.ch21.*tom.ch21);
    end
elseif nargin == 2
    S1 = varargin{1};
    S2 = varargin{2};
    dim = find(size(S1)==3);    
    intensity = sqrt(dot(S1,S1,dim)) + sqrt(dot(S2,S2,dim));
end




