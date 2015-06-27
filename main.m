clear
clc
close all
%% Apply EKF to random-set result from GrpLasso
% --- \sigma_e^2 = 0.265338 : Observation noise

%load hpcc_grpLasso_0724149am
load alien_0213
% --- Load calibration matrix
load calibmatrix calib
Ma=calib(1:2,1:2);
Mb=calib(3,1:2);
Mc=calib(1:2,3);
Mk=calib(3,3);

itmp = (1:2);
jtmp = kron((1:1)',[1;1]);

% --- Load the input command
[~, ~, raw] = xlsread('.\FeaturesAndLocations.xlsx','Sheet1','A2:AL464');
data = reshape([raw{:}],size(raw));
%----------------------------Initialize------------------------------------
dataSize = 400;
scaleRatio = 1;
input = downSampling([data(1:dataSize,1) data(1:dataSize,2)],scaleRatio);
time_dsamp=downSampling(data(1:dataSize,5),scaleRatio);
delta_t=zeros(size(time_dsamp,1),1);
mu_y = mean(data(1:dataSize,3:4));
var_y = sqrt(var(data(1:dataSize,3:4)));

for i=1:(size(time_dsamp,1)-1)
    delta_t(i)= time_dsamp(i+1)-time_dsamp(i);
end

dataSize = size(input,1);

[IC,IX] = sort(index(validateIndex+1:testIndex));
y_test_sort = y_test(IX',:);
y_est_test_sort = y_guess_test(IX',:);

% --- Create vector indicates observations of train, validate and test
% --- 0 :train, 1: validate, 2: test
obs_ix = zeros(dataSize,1);
for i=1:dataSize 
   if (isempty(find(IC==i,1))==0)
        obs_ix(i) = 1;
    end
end

phi = 0; % Heading angle
y_(1,1) = y(1,1);
y_(1,2) = y(1,2);

for i=2:dataSize
    %% Without EKF
    % --- Convert to world coord:
    y_pre(1,1) = y_(i-1,1)*var_y(1) + mu_y(1);
    y_pre(2,1) = y_(i-1,2)*var_y(2) + mu_y(2);
    y_w_temp = (Ma-y_pre*Mb)^(-1)*(y_pre*Mk-Mc);
    y_w_pre = y_w_temp';
    
    % --- Update position in world coord
    % --- input(i,1): linear speed
    % --- input(i,2): turn rate
    phi = phi + input(i-1,2)*delta_t(i-1);
    stmp = reshape([cos(phi),sin(phi)]',[],1);
    F = sparse(itmp,jtmp,stmp,2,1);
    F_ = full(F);
    delta_y_w = delta_t(i-1)*F_ * input(i-1,1);
    y_w  = y_w_pre + delta_y_w';
    
    % --- Convert back to frame coord
    y_f_temp = [y_w(1);y_w(2);1];
    fTmp = calib*y_f_temp;
    y_(i,1)=((fTmp(1)/fTmp(3))-mu_y(1))/var_y(1);
    y_(i,2)=((fTmp(2)/fTmp(3))-mu_y(2))/var_y(2);
    
end
clear y_pre phi delta_y_w y_w
% --- Initialize EKF
%model_uncertainty = 0.7039; % RMSE of y and y_
%observation_noise = 0.5151; % RMSE of Group LASSO (y_test and y_guess_test)
model_uncertainty = (0.02)^2;
observation_noise = (0.5151)^2;
Sig_w = [model_uncertainty 0 0;...
    0 model_uncertainty 0;...
    0 0 0.0];

P = eye(3)*0.01; % Initial states are EXACTLY same as true, so covariance matrix has small value
M = [1 0 0;
    0 1 0];
R = [observation_noise 0;
    0 observation_noise]; % measurement covariance

yhat = y(1,:);
hhat = 0;
y_ekf(1,1) = y(1,1);
y_ekf(1,2) = y(1,2);
obv_count = 1;
P_rec(:,:,1) = P;
figure(1)

square_error_ekf = 0;
for i=2:dataSize
    %% --- EKF prediction:
    % --- Convert to world coord
    y_pre(1,1) = y_ekf(i-1,1)*var_y(1) + mu_y(1);
    y_pre(2,1) = y_ekf(i-1,2)*var_y(2) + mu_y(2);
    y_w_temp = (Ma-y_pre*Mb)^(-1)*(y_pre*Mk-Mc);
    y_w_ekf_pre = y_w_temp';
    
    % --- Update position in world coord
    hhat_ = hhat + input(i-1,2)*delta_t(i-1);
    stmp = reshape([cos(hhat),sin(hhat)]',[],1);
    F = sparse(itmp,jtmp,stmp,2,1);
    F_ = full(F);
    delta_y_w = delta_t(i-1)*F_ * input(i-1,1);
    y_w  = y_w_ekf_pre + delta_y_w';
    
    % --- Convert back to frame coord
    y_f_temp = [y_w(1);y_w(2);1];
    fTmp = calib*y_f_temp;
    y_hat_(1,1)=((fTmp(1)/fTmp(3))-mu_y(1))/var_y(1);
    y_hat_(1,2)=((fTmp(2)/fTmp(3))-mu_y(2))/var_y(2);
    
    delta_F = [1 0 -input(i-1,1)*sin(hhat);...
        0 1  input(i-1,1)*cos(hhat);...
        0 0  1];
    P_ = delta_F*P*delta_F' + Sig_w;
    %% --- EKF Update
    % [K] = 3\times 2
    if (obs_ix(i)==1) % observation is made
        %ix_temp = find(IC==i);
        loc_noise = y_est_test_sort(obv_count,:);
        K = P_*M'*(M*P_*M'+R)^-1;
        temp_hat_ = [y_hat_ hhat_]';
        temp_hat = temp_hat_ + K*(loc_noise' - M*temp_hat_);
        P = (eye(3) - K*M)*P_;
    
        obv_count = obv_count + 1;

        y_ekf(i,:) = temp_hat(1:2);
        hhat = temp_hat(end);
        % --- Update the RMSE of grpLASSO+EKF for test ONLY
        square_error_ekf = square_error_ekf + sum((y(i,:)-y_ekf(i,:)).^2);
    else
        P = P_;
        y_ekf(i,:) = y_hat_;
        hhat = hhat_;
    end
    P_rec(:,:,i) = P;
    % --- Plot error ellipse
    %error_ellipse('mu',y(i,:),'C',P(1:2,1:2),'conf',0.9);
end
figure(1)
hold on
fprintf('RMSE of pure kinematics: %f\n',sqrt(mseCal(y,y_)));
fprintf('RMSE of LASSO: %f\n',sqrt(mseCal(y_est_test_sort,...
                                                            y_test_sort)));
fprintf('RMSE of LASSO+EKF (full): %f\n',sqrt(mseCal(y,y_ekf)));
fprintf('RMSE of LASSO+EKF (test only): %f\n',sqrt(square_error_ekf/(obv_count)));

plot(y(:,1),y(:,2),'r-','LineWidth',2.5,'MarkerSize',7);
plot(y_(:,1),y_(:,2),'k:','LineWidth',2.5,'MarkerSize',7);
plot(y_ekf(:,1),y_ekf(:,2),'b-.','LineWidth',2.5,'MarkerSize',7);
plot(y_est_test_sort(:,1),y_est_test_sort(:,2),'gs','LineWidth',2.5);
hold off

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);
xlabel('Y axis (meters)','FontSize',18);
ylabel('X axis (meters)','FontSize',18);
grid on
box on

figure(2)
hold on
line = 0:0.01:0.05;
for i=20:80
    if (i>1)
        plot([i-1 i],[P_rec(1,1,i-1) P_rec(1,1,i)],'k-','LineWidth',2);
        %plot(i,P_rec(1,2,i),'b*','MarkerSize',5);
        %plot(i,P_rec(2,1,i),'gd','MarkerSize',5);
        plot([i-1 i],[P_rec(2,2,i-1) P_rec(2,2,i)],'k-.','LineWidth',2);
        if (isempty(find(IC==i, 1))==0)
            plot(ones(size(line,2),1)*i,line,'k:');
        end
    end
end
xlim([19 80])
hold off;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);
xlabel('Iteration','FontSize',18);
ylabel('P_{[i,i]}','FontSize',18);
box on

% figure(3)
% lowBound = 1;
% upBound = 400;
% lineSpace = 1;
% scale = 1;
% 
% plot_index = lowBound:lineSpace:upBound;
% %plot_index = plot_index' + ones(size(plot_index,1),1)*obs_num;
% hold on
% plot(scale*y(plot_index,1),scale*y(plot_index,2),'kd-','LineWidth',2.5,'MarkerSize',7);
% plot(scale*y_ekf(plot_index,1),scale*y_ekf(plot_index,2),'kd-.','LineWidth',2.5,'MarkerSize',7);
% 
% for i = lowBound:lineSpace:upBound
%     temp = find(IC==i);
%     if (isempty(temp)==0)
%         plot(scale*y_est_test_sort(temp,1),scale*y_est_test_sort(temp,2),'ks','LineWidth',2.5);
%         plot(scale*y(i,1),scale*y(i,2),'rd','LineWidth',2.5,'MarkerSize',7);
%     end
%     error_ellipse('mu',scale*y(i,:),'C',P_rec(1:2,1:2,i),'conf',0.9);
% end
% hold off



