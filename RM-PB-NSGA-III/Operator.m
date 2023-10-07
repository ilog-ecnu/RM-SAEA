function OffspringDec = Operator(PopDec,M,Problem,K)
% Generate offsprings by the models

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------



    %% Parameter setting
    if nargin < 4
        K = 5;
    end
    [N,D]  = size(PopDec);
    %% Modeling
    [Model,probability,partition] = LocalPCA(PopDec,M,K);

    %% Reproduction
    OffspringDec = zeros(N,D);
    % Generate new trial solutions one by one
    for t = 1 : 50
         % Select one cluster by Roulette-wheel selection
        i = mod(t,N)+1;
        k = partition(i);
        MatingPool = find(partition==k);
        if length(MatingPool)<5
            MatingPool = 1:N;
        end
        % Generate one offspring
        if ~isempty(Model(k).eVector)
            lower = Model(k).a - 0.25*(Model(k).b-Model(k).a);
            upper = Model(k).b + 0.25*(Model(k).b-Model(k).a);
            trial = rand(1,M-1).*(upper-lower) + lower;
            Sample = Model(k).mean + trial*Model(k).eVector(:,1:M-1)' ;
            Q = randsample(MatingPool,2);
            OffspringDec(i,:) = MYDE(Problem,Sample,PopDec(Q(1),:),PopDec(Q(2),:));
        else
            OffspringDec(i,:) = Model(k).mean + randn(1,D);
        end
    end
    %Offspring = Problem.Evaluation(OffspringDec);
end