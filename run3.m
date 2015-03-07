%% Many fits to Amy's data
plate = growthFit5('plate1trial1.csv'); % generate a cell array with the data

% Generate a CSV table with all the data
fileID = fopen('plateOutput.csv','w');
fprintf(fileID,'well number,strain,logistic R^2,logistic doubling time [min],logistic max growth rate [per hr],logistic carrying capacity,logistic growth infleciton time [hr],gompertz R^2,gompertz max growth rate,gompertz carrying capacity,gompertz lag time [hr]\r\n');
for i=1:96
    fprintf(fileID,'%d,%s,%s,%2.7f,%2.7f,%2.7f,%2.7f,%2.7f,%2.7f,%2.7f,%2.7f,%2.7f',i,plate{i}.mutant,...
        plate{i}.logistic.R2,plate{i}.logistic.doublingTime,plate{i}.logistic.maxGrowthRate,plate{i}.logistic.carryingCapacity,plate{i}.logistic.growthInflection,plate{i}.gompertz.R2,plate{i}.gompertz.maxGrowthRate,plate{i}.gompertz.carryingCapacity,plate{i}.gompertz.lagTime);
    fprintf(fileID,'\r\n');
end
fclose(fileID);