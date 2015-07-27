% --- Plot trajectory
plot(y_guess_test_sort(:,1),y_guess_test_sort(:,2),'b*--','LineWidth',2);
hold on
plot(y_test_sort(:,1),y_test_sort(:,2),'rd--','LineWidth',2);

% --- Plot axies INDIVIDUALLY
subplot(2,1,1)
hold on
plot(y_guess_test(:,1),'b*--','LineWidth',2);
plot(y_test(:,1),'rd--','LineWidth',2);
title('X axis')
hold off
subplot(2,1,2)
hold on
plot(y_guess_test(:,2),'b*--','LineWidth',2);
plot(y_test(:,2),'rd--','LineWidth',2);
title('Y axis')
hold off

% --- Plot evolution of vector b
load resultGIST022011h00am.mat
figure(1)
hold on
for i=1:size(Bx,1)
	plot(FitInfoX.Lambda,Bx(i,:),'r','LineWidth',2);
	%plot(Bx(i,:),'r','LineWidth',2);
end
hold off
grid on
box on
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);
xlabel('Lambda','fontsize',18);
ylabel('B','fontsize',18);

figure(2)
plot(errorOfXTest,'LineWidth',2)
box on
grid on
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);
xlabel('Lambda','fontsize',18);
ylabel('RMSE','fontsize',18);


% --- plot trajectory 2 sets: train, validate and test
hold on
plot(y_est_train_sort(:,1),y_est_train_sort(:,2),'b*--','LineWidth',2);
plot(y_est_validate_sort(:,1),y_est_validate_sort(:,2),'b*--','LineWidth',2);
plot(y_est_test_sort(:,1),y_est_test_sort(:,2),'b*--','LineWidth',2);

%plot(y_guess_train{1}(:,1),y_guess_train{1}(:,2),'b*','LineWidth',2);
%plot(y_guess_validate{minIndex}(:,1),y_guess_validate{minIndex}(:,2),'bo','LineWidth',2);
%plot(y_guess_test(:,1),y_guess_test(:,2),'bs','LineWidth',2);

plot(y(:,1),y(:,2),'rd--','LineWidth',2);

% --- Plot coefficients of group LASSO

cm = colormap;
col_code = round(linspace(1,60,nf));
for i=1:size(B_rec,1)
	for j=1:size(B_rec,2)
		if(B_rec(i,j)<1e-3)
			B_rec(i,j) = 1e-3;
		end
	end
end
hold on
for i=1:3:nf
	plot(lambda(2:end),B_rec(i,:),'Color',cm(col_code(i),:),'LineWidth',1.5);
	plot(lambda(2:end),B_rec(i+30,:),'Color',cm(col_code(i),:),'LineWidth',1.5);
end
xlim([0.003 0.019]);
xlabel('Lambda','FontSize',14);
ylabel('Value of entries of matrix B','FontSize',14);
set(gca,'FontSize',16);
box on

% --------------- Plot survival features
for i =1:size(BOpt,1)
	if (BOpt(i,1)<1e-2)
		BOpt(i,1) = 0;
	end
end
bar(BOpt);
xlim([0 200]);
xlabel('Entries of the optimal matrix B','FontSize',14);
ylabel('Values','FontSize',14);
set(gca,'FontSize',16);
fprintf('Number of survived features: %d',(nnz(BOpt)));

% ------ Plot y_train and y_test
load hpcc_run_0403(2)_OP_SURF
hold on
plot(y_train(:,1),y_train(:,2),'k:','LineWidth',2);
plot(y_test(:,1),y_test(:,2),'k','LineWidth',2);
hold off
xlabel('(meters)','FontSize',14);
ylabel('(meters)','FontSize',14);
set(gca,'FontSize',16);
box on

