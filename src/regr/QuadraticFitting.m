%% QuadraticFitting.m
% fit data from CSV file 
% file format is as follows
% Y,X1,X2 ...
% The first line will be skipped as header 

%% 
function QuadraticFitting(filename)

% reading data 
% start from second line and first column 
M = csvread(filename, 1, 0);
y = M(:,1);
X = M(:,2:end);

% fitting data 
whichstats = {'beta', 'yhat', 'r'};
stats = regstats(y, X, 'quadratic', whichstats);

% output data 
% stats have different dimensions of report values 
% dump out them separately
[pathstr, name, ext] = fileparts(filename);
% write beta 
outfilename = fullfile(pathstr, [name '-beta' ext]);
T = table(stats.beta);
T.Properties.VariableNames = {'beta'};
writetable(T, outfilename, 'Delimiter', ',');
% write yhat and r
outfilename = fullfile(pathstr, [name '-y' ext]);
T = table(stats.yhat, stats.r);
T.Properties.VariableNames = {'yhat', 'r'};
writetable(T, outfilename, 'Delimiter', ',');
end