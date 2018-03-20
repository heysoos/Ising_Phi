function [ J ] = make_J(N,Params,options)
% Make a connectivity matrix J.
% 
% Params: [sparsity, anti-correlated sparsity, seed]
% options: 'norm', 'normpos', 'ones':
%   'norm':     returns a normally distributed connectivity
%   'lognorm':  lognormally distributed connectivity
%   'normpos':  same as norm but values are between 0 and 1
%   'ones':     returns a matrix of ones .

if isnumeric(Params)
    spars = Params(1); % 0 means full connected, 1 means no connections
    antiCorSpars = 1;
    if size(Params,2) > 1
        antiCorSpars = Params(2); % 0 means all anticorrelations, 1 means no anticorrelations
    end
    if size(Params, 2) > 2
        if seed ~= 0
            seed = Params(3);
        end
        rng(seed)
        fprintf(['Setting seed to: ', num2str(seed), '.\n'])
    end
else
    error('Input numeric sparsity!')
end


sparsMap = (rand(N) > spars);
sparsMap = sparsMap.*((rand(N) > 1- antiCorSpars)*2 -1);
sparsMap(1:N+1:end) = 0;
sparsMap = triu(sparsMap) + triu(sparsMap,1)';

if ~(spars == 1)
    while (sum(sum(sparsMap)) == 0)
        sparsMap = (rand(N) > spars);
        sparsMap = sparsMap.*((rand(N) > 1- antiCorSpars)*2 -1);
        sparsMap(1:N+1:end) = 0;
    end
end

% if nargin == 3
%     if strcmp(options,'norm')
%         J = randn(N);
%         J = (J - min(min(J)))./(max(max(J)) - min(min(J)));
%         J = J.*sparsMap;
%     elseif strcmp(options,'ones')
%         J = ones(N);
%         J(1:N+1:end) = 0;
%     else
%         error('Option does not exist!')
%     end
% end

if nargin == 3
    switch options
        
        case 'norm'
            
            J = randn(N);
            J = J.*sparsMap;
            
        case 'normpos'
            
            J = randn(N);
            J = (J - min(min(J)))./(max(max(J)) - min(min(J)));
            J = J.*sparsMap;
            
        case 'lognorm'
            
            J = lognrnd(0,1, N, N);
            J = J./max(max(J));
            J = J.*sparsMap;
           
            
        case 'ones'
            
            J = ones(N);
            J(1:N+1:end) = 0;
            
        otherwise
            error('Invalid option. Try: ''norm'', ''normpos'', or ''ones''.')
    end
else
    J = rand(N);
    J = J.*sparsMap;
end

J = triu(J) + triu(J,1)';

end

