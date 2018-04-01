function [num_autoencoders,Abundance,Endmembers] = SNSAnew(Y,m,max_num_autoencoders,FmaxIters)
%%    Initialization
v= 'off';
decayP = 0;     %decay factor for positive weights 
decayN = 1;     %decay factor for negative weights 
lrateRO = 0.01;    
lrateIP = 0.001;   
regRO = 0.0002; 
meanIP = 0.02;
stop_threshold = 0;
Ynew=Y;
num_autoencoders =0;
%% first group of autoencoders 
while num_autoencoders <= max_num_autoencoders && stop_threshold == 0
    disp(['number of autoencoders' num2str(num_autoencoders) '/' num2str(max_num_autoencoders)]);% exhibit Iteration
    [EndmCandidates, indice, ~ ]= VCA(Ynew,'Endmembers',30,'verbose',v);
    [samples,~] = trainSample(Ynew,20,EndmCandidates,1); 
    N = size(samples,2);
    inpDim = size(samples{1},1); 
    hidDim = size(samples{1},1); 
for j1 = 1:N
  W = 0.5 * (rand(inpDim,hidDim) - 0.5 * ones(inpDim,hidDim));
  c = ones(hidDim,1);           
  d = -3 * ones(hidDim,1);
  numSamples = size(samples{j1},2);
  for e = 1:FmaxIters
         p = randperm(numSamples); % for randomized presentation 
      for i1=1:numSamples 
          inp = samples{j1}(:,p(i1));
          g = W' * inp;
          h = 1 ./ (1 + exp(-c .* g - d));
          out = W * h;
          lrate = lrateRO/(regRO + sum(h.^2)); 
          error = inp - out;
          W = abs(W) + lrate * error * h';
          if  decayP > 0
              idx = find(W > 0);
              W(idx) = W(idx) - decayP * W(idx);
          end                
          % decay function for negative weights
          if  decayN == 1
          % pure NN weights!
              W = max(W, 0);
          elseif decayN > 0
              idx = find(W < 0);
              W(idx) = W(idx) - decayN * W(idx);
          end                
             hones = ones(hidDim,1);
             tmp = lrateIP * (hones - (2.0 + 1.0/meanIP) * h + h.^2/meanIP);    
             c = c + tmp;
             d = d + lrateIP * hones ./ d + g .* tmp; 
      end 
      Xout{j1}=out;
  end
end
for k_1=1:N
    a1 = Xout{k_1};b1 = EndmCandidates(:,k_1);
    errRadians(k_1) = acos(dot(a1, b1)/ (norm(b1) * norm(a1)));
end
threshold1 = mean(errRadians)+3*std(errRadians,0,2);
threshold = min(threshold1,1);
minimi=0.04;
if threshold < minimi;
    stop_threshold = 1;
end
OutlierIndex=find(errRadians>=threshold);
outlier_index=indice(:,OutlierIndex);
Ynew(:,outlier_index)=[];
maxIndex=find(threshold1>=minimi);
lg=length(maxIndex);
for l=1:lg
    Ynew(:,1)=[];
end
num_autoencoders = num_autoencoders +1;
end
%% second group of autoencoders
disp('**** SNSA for unmixing  ****')
v='off';
[~,Pixels] = size(Ynew);
[InitialW, ~, ~ ]= VCA(Ynew,'Endmembers',m,'verbose',v);
InitialH = zeros(length(InitialW(1,:)),Pixels);
endm = [1e-5*InitialW;ones(1,length(InitialW(1,:)))];
for i2=1:Pixels
    Rmixed = [1e-5*Ynew(:,i2); 1];
    InitialH(:,i2) = lsqnonneg(endm,Rmixed);
end
InitialW=abs(InitialW);
W = InitialW;
H = InitialH;
[PCA, ~] = princomp(Ynew',0); 
RMSE = sum(sum((Ynew-W*H).^2)); 
RMSE = [RMSE, 0];
hint1 = 0; hint2 = 0;
for Iteration  = 1:500
    if RMSE(1)-RMSE(2)>0.0001 
       hint1 = 0;
    else
       hint1 = hint1+1;hint2 = hint2+1;
    end
    if Iteration < 5
       hint1 = 0;
    end
%  Assure abundance matrix H to Sum-to-One Constraint 
   R_augment = [Ynew; 20*ones(1,Pixels)];
   W_augment = [W; 20*ones(1,m)];
   Threshold = 0.0001;
   MaxIteration = 200;
   [H] = SNSAabun(R_augment, W_augment, H, MaxIteration, Threshold);
   D = PCA(:,1:m-1);
   lam = 0.014;
   
   RR=Ynew';
   meanData = mean(RR); 
   inp = Ynew; 
   PP = H*H';
   nn = -0.001;
   B = [ones(1,m); zeros(m-1,m)];
   C = [zeros(1,m-1); eye(m-1)]*D';
   meanData = meanData'*ones(1,m);
   I = B+C*(W-meanData);
   detz2 = det(I)^2;
   ID = pinv(I)*C;
   tWp = H*H'*W' - H*inp' + lam * detz2 * ID;
   conjp = tWp;
   W = W + nn*tWp';
    for iteration = 1:50
     I = B+C*(W-meanData);
     detz2 = det(I)^2;
     ID = pinv(I)*C;
     tWpp = H*H'*W' - H*inp' + lam * detz2 * ID;
     beta = abs(sum(tWpp.*(tWpp-tWp),1)./(sum(tWp.^2,1)));
     conj = -tWpp+ones(m,1)*beta.*conjp;
     AAd = H*H'*conj;
     alpha = sum(conj.*(-tWpp),1)./max(sum(conj.*AAd,1),eps);
     tW = W' + conj.*repmat(alpha,size(PP,2),1);
     W=tW';
 % decay function for negative weights of NNSAE
     if  decayP > 0
         idx = find(W > 0);
         W(idx) = W(idx) - decayP * W(idx);
     end               
     if  decayN == 1
         W = max(W, 0);
     elseif decayN > 0
         idx = find(net.W < 0);
         W(idx) = W(idx) - decayN * W(idx);
     end
     tWp = tWpp;
     conjp = conj;
    end
   ObNew = sum(sum((Ynew-W*H).^2));
   RMSE = [RMSE, ObNew];  
   if hint1>5 && hint2>20 
      break;
   end
end
Endmembers=W;% obtain endmembers
Abundance=H; % obtain abundance
%% selec training samples
function [samples,index] = trainSample(Y,k,EndmembersOriginal,test) 
% input __________________________________________________________________
% Y is hyperspectral data
% k is the number of per class sample for training
% p is the number of endmembers
% output __________________________________________________________________
% samples is traing samples
% By Yuanchao Su, 2016.

[~, Nt] = size(Y);
U = EndmembersOriginal;    
switch test
case 1
     for i=1:Nt
         for j=1:size(U,2)
              Ht(j,i)=norm((U(:,j)-Y(:,i)),2);
         end
     end
Bc = sort(Ht,2);
     for i= 1:size(Bc,1) 
         Z=Bc(i,:);     
         Q=Z(1:k); 
         index(i,:)=find(Ht(i,:)<=Q(k));
         samples{1,i}=Y(:,index(i,:));
     end
[mt,nt]=size(index);
index=reshape(index',1,mt*nt);

case 2  
    for i= 1:Nt
        for j=1:size(U,2)
            a = U(:,j); b = Y(:,i);
            Ht(j,i) = acos(dot(a, b)/ (norm(b) * norm(a)));
        end
    end
Bc = sort(Ht,2);
     for i= 1:size(Bc,1) 
         Z=Bc(i,:);     
         Q=Z(1:k); 
         index(i,:)=find(Ht(i,:)<=Q(k));
         samples{1,i}=Y(:,index(i,:));
     end
[mt,nt]=size(index);
index=reshape(index',1,mt*nt);
samples = cell2mat(samples);

case 3
Ymean = mean(Y');
Ymean = Ymean';
     for i=1:Nt
         Ht(i)=norm((Ymean-Y(:,i)),2);% Euchliod distance
     end
     Bc = sort(Ht,2);
     Q = Bc(1:k);
     index = find(Ht <=Q(k));
     samples = Y(:,index);
end
end
end
