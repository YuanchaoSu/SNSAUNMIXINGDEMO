clear;
clc;
close all
%% Test it with different synthetic data
% There are 4 hyperspectral data for this test
% set test=1, 2, 3 or 4
test = 1;

%% Data 1
switch test
case 1
[Y, X, M]=LinearSyn(1);
FmaxIters=500;
max_num_autoencoders=15;
m=4;
[num_autoencoders,~,Endmembers] = SNSAnew(Y,m,max_num_autoencoders,FmaxIters);

% test N-FINDR
[bands,np]=size(Y);
[numBands,Newpixels]=size(Y);
img = reshape(Y.', 1, Newpixels, numBands); 
[nfindrEndms,~] = NFINDR(img,m);
% test VCA
v='off';
[A_vca, ~, ~ ]= VCA(Y,'Endmembers',m,'verbose',v);
vcaEndms=A_vca;
% display
figure
set(gcf, 'position', [400 400 1800 500]);

subplot(1,3,1);
hold on
plot(nfindrEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('NFINDR');
hold off

subplot(1,3,2);
hold on
plot(vcaEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('VCA');
hold off

subplot(1,3,3);
hold on
plot(Endmembers,'r','LineWidth',2);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',2);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('SNSA');
hold off
%% Data 2
case 2
[Y, X, M]=LinearSyn(2);
FmaxIters=500;
max_num_autoencoders=15;
m=4;
[num_autoencoders,~,Endmembers] = SNSAnew(Y,m,max_num_autoencoders,FmaxIters);

% test N-FINDR
[bands,np]=size(Y);
[numBands,Newpixels]=size(Y);
img = reshape(Y.', 1, Newpixels, numBands); 
[nfindrEndms,~] = NFINDR(img,m);
% test VCA
v='off';
[A_vca, ~, ~ ]= VCA(Y,'Endmembers',m,'verbose',v);
vcaEndms=A_vca;
% display
figure
set(gcf, 'position', [400 400 1800 500]);

subplot(1,3,1);
hold on
plot(nfindrEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('NFINDR');
hold off

subplot(1,3,2);
hold on
plot(vcaEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('VCA');
hold off

subplot(1,3,3);
hold on
plot(Endmembers,'r','LineWidth',2);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',2);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('SNSA');
hold off
%% Data 3
case 3
[Y, X, M]=LinearSyn(test);
FmaxIters=500;
max_num_autoencoders=15;
m=4;
[num_autoencoders,Abundance,Endmembers] = SNSAnew(Y,m,max_num_autoencoders,FmaxIters);
% test N-FINDR
[bands,np]=size(Y);
[numBands,Newpixels]=size(Y);
img = reshape(Y.', 1, Newpixels, numBands); 
[nfindrEndms,~] = NFINDR(img,m);
% test VCA
v='off';
[A_vca, ~, ~ ]= VCA(Y,'Endmembers',m,'verbose',v);
vcaEndms=A_vca;
% display
figure
set(gcf, 'position', [400 400 1800 500]);

subplot(1,3,1);
hold on
plot(nfindrEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('NFINDR');
hold off

subplot(1,3,2);
hold on
plot(vcaEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('VCA');
hold off

subplot(1,3,3);
hold on
plot(Endmembers,'r','LineWidth',2);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',2);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('SNSA');
hold off
%%
case 4
    
[Y, X, M]=LinearSyn(test);
FmaxIters=500;
max_num_autoencoders=15;
m=4;
[num_autoencoders,Abundance,Endmembers] = SNSAnew(Y,m,max_num_autoencoders,FmaxIters);
% test N-FINDR
[bands,np]=size(Y);
[numBands,Newpixels]=size(Y);
img = reshape(Y.', 1, Newpixels, numBands); 
[nfindrEndms,~] = NFINDR(img,m);
% test VCA
v='off';
[A_vca, ~, ~ ]= VCA(Y,'Endmembers',m,'verbose',v);
vcaEndms=A_vca;
% display
figure
set(gcf, 'position', [400 400 1800 500]);

subplot(1,3,1);
hold on
plot(nfindrEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('NFINDR');
hold off

subplot(1,3,2);
hold on
plot(vcaEndms,'r','LineWidth',3);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',3);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('VCA');
hold off

subplot(1,3,3);
hold on
plot(Endmembers,'r','LineWidth',2);
plot(M,'-.','color',[0,0.45,0.74],'LineWidth',2);
ylabel('reflectance','fontsize',20);
xlabel('bands','fontsize',20);
set(gca,'FontSize',15);
title('SNSA');
hold off
end