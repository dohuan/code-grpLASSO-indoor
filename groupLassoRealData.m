clear
clc
close all
%%
tic
addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))
[~, ~, raw] = xlsread('.\FeaturesAndLocations.xlsx','Sheet1','A2:AL464');
data = reshape([raw{:}],size(raw));
%----------------------------Initialize------------------------------------
dataSize=400;
scaleRatio = 1;
nf=30;  %number of features, maximum: 32
dim=2; %number of member of each group
iter=30; % iteration times for lambda
%index=randperm(6);
b=sdpvar(nf,1);
sdpvar t;
B=sdpvar(dim*nf,1);
error=zeros(iter,2);
lambdaMin=0;
lambdaStep=0.0001;
%-----------------------------Load data------------------------------------
y=downSampling([data(1:dataSize,3) data(1:dataSize,4)],scaleRatio);
f=[];
for i=1:nf
    f=[f,downSampling(data(1:dataSize,i+6),scaleRatio)]; %data starts at 7th column
end
%f=[data(1:dataSize,7) data(1:dataSize,8) data(1:dataSize,9)...
%    data(1:dataSize,10) data(1:dataSize,11) data(1:dataSize,12)];

%-------------------------normalize features and targets-------------------
for i=1:nf
   f(:,i)=zscore(f(:,i)); 
end
y(:,1) = zscore(y(:,1));
y(:,2) = zscore(y(:,2));
nt=size(y,1); % post-downSampled data size is nt
lambda = lambdaMin:lambdaStep:(lambdaMin+lambdaStep*iter);
%lambda=zeros(iter);
B_rec=zeros(nf*2,iter);
I=eye(dim,dim);
% --- Update dataSize
dataSize = size(y,1);

%----divide data into TRAIN, VALIDATE and TEST sets------------------------

trainIndex = round(dataSize/2);
validateIndex = round(dataSize*3/4);
testIndex = dataSize;
mode = 0; % 1 for ordered, 0 for random
if (mode==1)
    y_train = y(1:trainIndex,:);
    f_train = f(1:trainIndex,:);
    y_validate = y(trainIndex+1:validateIndex,:);
    f_validate = f(trainIndex+1:validateIndex,:);
    y_test = y(validateIndex+1:testIndex,:);
    f_test = f(validateIndex+1:testIndex,:);
else
    index = randperm(nt);
    y_train = y(index(1:trainIndex),:);
    f_train = f(index(1:trainIndex),:);
    y_validate = y(index(trainIndex+1:validateIndex),:);
    f_validate = f(index(trainIndex+1:validateIndex),:);
    y_test = y(index(validateIndex+1:testIndex),:);
    f_test = f(index(validateIndex+1:testIndex),:);
end
%------------- Train the group LASSO with TRAIN set -----------------------
yhat=reshape(y_train,[],1);
xhat=kron(I,f_train); % (I \kronecker x)\times B

F=[];
for i=1:nf
    W_temp = sparse([1;2],[i;i+nf],ones(2,1),2,2*nf);
    F=[F,cone(W_temp*B,b(i))];
end
t=sum(b);
F=[F,t,set(B>0)];

for i=1:iter
    fprintf('process: %d/%d',i,iter);
    solvesdp(F,1/(2*dataSize)*norm(yhat-xhat*B,2)+lambda(i)*t);
    y_guess_train{i} = f_train*reshape(value(B),[],2);
    error_train(i,1) = norm(y_train(:,1)-y_guess_train{i}(:,1));
    error_train(i,2) = norm(y_train(:,2)-y_guess_train{i}(:,2));
    B_rec(:,i)=value(B);    
end

stopWatch=toc;
s=sprintf('\n time ellapse for training: %f minutes',stopWatch/60);
disp(s);
%------------- Use VALIDATE set to find optimal lambda --------------------
valSize = size(y_validate,1);
y_guess_validate = cell(iter,1);
for i=1:iter
    y_guess_validate{i,1} = f_validate*reshape(B_rec(:,i),[],2);
    error_validate(i,1) = mse(y_guess_validate{i}-y_validate);
end
minIndex = find(error_validate == min(error_validate));
BOpt = B_rec(:,minIndex);
%--------------- Apply to TEST set ----------------------------------------

y_guess_test = f_test*reshape(BOpt,[],2);


figure('name','error train');
subplot(2,1,1);
plot(error_train(:,1));
subplot(2,1,2);
plot(error_train(:,2));

figure('name','error validate');
plot(error_validate,'LineWidth',2);

figure('name','test X');
plot(y_test(:,1),'LineWidth',2);
hold on;
plot(y_guess_test(:,1),'Color','red','LineWidth',2);
%plot(y(:,1),'Color','green');
hold off;

figure('name','test Y');
plot(y_test(:,2),'LineWidth',2);
hold on;
plot(y_guess_test(:,2),'Color','red','LineWidth',2);
%plot(y(:,2),'Color','green');
hold off;