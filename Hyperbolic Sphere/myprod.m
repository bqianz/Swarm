function res = myprod(varargin)
if nargin > 1
    S = size(varargin{1});
    res = ones(S);
    for i=1:nargin
        if varargin{i} == zeros(S)
            res = zeros(S);
            break;
        else
            res = res.*varargin{i};
        end
    end
else
    disp('Not enough inputs.');
end
end