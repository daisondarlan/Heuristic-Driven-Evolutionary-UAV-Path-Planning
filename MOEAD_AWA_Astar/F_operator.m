function newpop = F_operator(population,MatingPool,Boundary,model,kf)
% This function generates a new population by genetic operators
    
    N = size(MatingPool,2);
    D = 3;
%-----------------------------------------------------------------------------------------
% Parameters setting
    ProC = 1;       % The probability of crossover
    ProM = 1/D;     % The probability of mutation
    DisC = 20;   	% The parameter of crossover
    DisM = 30*kf;   	% The parameter of mutation
%-----------------------------------------------------------------------------------------
% Simulated binary crossover
    Parent1 = [];
    Parent2 = [];
    Parent1 = [MatingPool(1).rnvec];
    Parent2 = [MatingPool(2).rnvec];
    wp = randi([2,size(Parent1,1)-1]);

    ptt = aStarCrossover(Parent1,Parent2,model);
    if isempty(ptt)
        Offspring = [Parent1(1:wp,:);Parent2(wp+1:end,:)];
    else
        Offspring = ptt;
    end

%-----------------------------------------------------------------------------------------
% Polynomial mutation
    if rand<1 %Using the DTLZ mutation strategy
        MaxValue = repmat(Boundary(1,:),N*size(population(1).rnvec,1),1);
        MinValue = repmat(Boundary(2,:),N*size(population(1).rnvec,1),1);
        k    = rand(size(population(1).rnvec,1),D);
        miu  = rand(size(population(1).rnvec,1),D);
        Temp = k<=ProM & miu<0.5;
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
        Temp = k<=ProM & miu>=0.5; 
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
    %-----------------------------------------------------------------------------------------
        rnvec_mut = aStarMutation(Offspring,model);
        Offspring = rnvec_mut;
        newpop = population;
        newpop.rnvec = Offspring;
        newpop.rnvec = sortrows(newpop.rnvec,randi(2));
        newpop.path = testBspline([newpop.rnvec(:,1)';newpop.rnvec(:,2)'],model.xmax)';            
        newpop = adjust_constraint_turning_angle(newpop,model);            
        newpop = evaluate(newpop);
    else %Using no mutation strategy, just crossover
        newpop = population;
        
    end


end