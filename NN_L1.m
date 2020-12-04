clc;clear all;close all
%%
warning('off','all');
rmpath(genpath(sprintf('%s%s',pwd)));
addpath(genpath(sprintf('%s%s',pwd,'/L1GeneralExamples')))
warning('on','all');
% addpath(sprintf('%s%s',pwd,loc));
%%
fig = 100;
lambda = 1;
options.quadraticInit = 1; % Use quadratic initialization of line search
options.verbose = 1;
% 
% % Generate non-linear regression data set
% nInstances = 200;
% nVars = 1;
% [X,y] = makeData('regressionNonlinear',nInstances,nVars);
load 'training_data2.mat'
max_ati = max(abs(atiTrain));

mean_x = mean(xTrain);
stddev_x = std(xTrain);
xmn = meanNormalize(xTrain,mean_x,stddev_x);
atin = atiTrain./max_ati;

xmnt = xmn(1:floor(0.9*size(xTrain,1))-1,:);
atit = atin(1:floor(0.9*size(atiTrain,1))-1,:);

xmnValid = xmn(floor(0.9*size(xTrain,1)):end,:);
atiValid = atin(floor(0.9*size(atiTrain,1)):end,:);

X = xmnt(1:50:end,:);
y = atit(1:50:end,:);

Xtest = xmnValid;
ytest = atiValid;

nInstances = size(X,1);
nVars = size(X,2);

X = [ones(nInstances,1) X];
nVars = nVars+1;

% Train neural network w/ multiple hiden layers
nHidden = [5 5];
nParams = nVars*nHidden(1);
for h = 2:length(nHidden)
    nParams = nParams+nHidden(h-1)*nHidden(h);
end
nParams = nParams+nHidden(end)*size(y,2);

funObj = @(weights)MLPregressionLoss_efficient(weights,X,y,nHidden);
fprintf('Training neural network for regression...\n');
lambdaL2 = 1e-3;
wMLP = randn(nParams,1);
funObjL2 = @(w)penalizedL2(w,funObj,lambdaL2);
while 1
    w_old = wMLP;
    wMLP = L1General2_PSSas(funObjL2,wMLP,lambda*ones(nParams,1),options);
    if norm(w_old-wMLP,inf) < 1e-5
        break;
	end
end
Xtest = [ones(size(Xtest,1),1) Xtest];
%% Plot results
figure;
hold on

yhat = MLPregressionPredict_efficient(wMLP,Xtest,nHidden,size(ytest,2));
plotLocal(ytest,yhat)
% h=plot(yhat,'g-');
% set(h,'LineWidth',3);
% legend({'Data','Neural Net'});
%%
% Form weights
inputWeights = reshape(wMLP(1:nVars*nHidden(1)),nVars,nHidden(1));
offset = nVars*nHidden(1);
for h = 2:length(nHidden)
    hiddenWeights{h-1} = reshape(wMLP(offset+1:offset+nHidden(h-1)*nHidden(h)),nHidden(h-1),nHidden(h));
    offset = offset+nHidden(h-1)*nHidden(h);
end
outputWeights = wMLP(offset+1:offset+nHidden(end));

% Make adjacency matrix
adj = zeros(nVars+sum(nHidden)+1);
for i = 1:nVars
    for j = 1:nHidden(1)
        if abs(inputWeights(i,j)) > 1e-4
            adj(i,nVars+j) = 1;
        end
    end
end
for h = 1:length(nHidden)-1
    for i = 1:nHidden(h)
        for j = 1:nHidden(h+1)
            if abs(hiddenWeights{h}(i,j)) > 1e-4
                adj(nVars+sum(nHidden(1:h-1))+i,nVars+sum(nHidden(1:h))+j) = 1;
            end
        end
    end
end
for i = 1:nHidden(end)
    if abs(outputWeights(i)) > 1e-4
        adj(nVars+sum(nHidden(1:end-1))+i,end) = 1;
    end
end

labels = cell(length(adj),1);
for i = 1:nVars
    labels{i,1} = sprintf('x_%d',i);
end
for h = 1:length(nHidden)
    i = i + 1;
    labels{i,1} = sprintf('b_%d',h);
    for j = 2:nHidden(h)
        i = i + 1;
        labels{i,1} = sprintf('h_%d_%d',h,j);
    end
end
labels{end,1} = 'y';


function xm = meanNormalize(x,mean_x,stddev_x)
xm = (x-mean_x)./stddev_x;
end


function plotLocal(x,y)
index = [1,3,5,2,4,6];
for i=1:6
    subplot(3,2,index(i))
    plot(x(:,i));hold on; plot(y(:,i))
end
end