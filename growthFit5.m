% Gompertz + Logistic Growth curve fitting
% Zack Scholl 2011-2013  zns@duke.edu
% Models are based of ZWIETERING 1990 "Modeling of the Bacterial Growth
% Curve"
% Takes input file formated as follows:
% CSV: First column is time in rows 3 - end
% All other columns are data points in row 3 - end and row 1+2 are names
function [plate] =  growthFit5(file)

c = importdata(file);
strainNames = c.textdata{1};
tmp = regexp(strainNames,'([^ ,:]*)','tokens');
strainNames = cat(2,tmp{:});
mutantNames = c.colheaders;
mutantNames(:,1)=[];time = c.data(:,1);
wells = c.data;
wells(:,1) = [];
T=time;
T=T-T(1);


for wellNumber =1:size(wells,2)

    M = wells(:,wellNumber);
    M=M-min(M);
    plate{wellNumber}.time = T;
    plate{wellNumber}.data = M;
    plate{wellNumber}.mutant = mutantNames{wellNumber};
    

t=linspace(T(1),T(end)); 

% Logistic
    minSearch=fminsearch(@logisticError,[0.1,.1],[],T,M);
    r=minSearch(1);
    t0=minSearch(2);
    H=1./(1+exp(-r*(T-t0)));
    K=(H'*M)/(H'*H);
    yLogistic = K./(1+exp(-r*(t-t0)));
    yLogisticResidual = K./(1+exp(-r*(T-t0)))-M;
    fit = real(yLogistic)';
    sse = sum((real(yLogisticResidual)).^2);
    sst = sum((real(M-mean(M))).^2);
    R2 = 1 - abs(sse/sst); 
    um = r*K/4;
    tdbl = log(2)/(4*um/K);
    plate{wellNumber}.logistic.R2 = R2;
    plate{wellNumber}.logistic.curveFit = [t' fit];
    plate{wellNumber}.logistic.doublingTime = tdbl*60;
    plate{wellNumber}.logistic.maxGrowthRate = um;
    plate{wellNumber}.logistic.carryingCapacity = K;
    plate{wellNumber}.logistic.growthInflection = t0;
    
  
    % Gompertz
    options=optimset('fminsearch');
    options= optimset(options,'MaxFunEvals',1e99);
    minSearch=fminsearch(@gompertzError,[.3,1],options,T,M);
    b=minSearch(1);
    c=minSearch(2);
    H=exp(-exp(minSearch(1)-minSearch(2).*T));
    K=(H'*M)/(H'*H);
    um = K*c/exp(1);
    lambda = (b-1)/c;
    tdbl = log(2)/(4*um/K); %I'm not 100% sure about this one....
    yGompertz = K*exp(-exp(b-c.*t));
    yGompertzResidual =  K*exp(-exp(b-c.*T))-M;
    fit = real(yGompertz)';
    sse = sum((real(yGompertzResidual)).^2);
    sst = sum((real(M-mean(M))).^2);
    R2 = 1 - abs(sse/sst);
    plate{wellNumber}.gompertz.R2 = R2;
    plate{wellNumber}.gompertz.carryingCapacity = K;
    plate{wellNumber}.gompertz.lagTime = lambda;
    plate{wellNumber}.gompertz.maxGrowthRate = um;
    plate{wellNumber}.gompertz.curveFit = [t' fit];
    
    
end



end
