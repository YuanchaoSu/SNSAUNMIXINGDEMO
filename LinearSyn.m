function [Y, X, M]=LinearSyn(test);
% This function can produce a sythetic hyperspectral datum
switch test
case 1 
load matlabM
M = True1_M;
case 2
load matlabM2
M = True2_M;
case 3
load matlabM3
M = True3_M;
case 4
load matlabM4
M = True4_M;
end
%%
bands=224;
MaxPurity=0.8;
dimension = 64;
p=4;
label = ones((dimension/8)^2,1);
num = floor(length(label)/p);
  for i=1:p-1
      label((i-1)*num+1:i*num) = (i+1); 
  end
  ridx = randperm(length(label));
  label = label(ridx)';
  label = reshape(label,dimension/8,dimension/8);
  X = zeros(dimension,dimension,p);
  img = zeros(dimension,dimension);
  for i=1:dimension
      for j=1:dimension
          for cls = 1:p
              if label(floor((i-1)/8)+1,floor((j-1)/8)+1) == cls
                 tmp = zeros(p,1);
                 tmp(cls) = 1;
                 X(i,j,:) = tmp;
                 img(i,j) = p;
              end
          end
      end 
  end
  win = 7;
  H = ones(win,win)/(win*win);
  for i=1:p
      X(:,:,i) = filter2(H,X(:,:,i));
  end
  X = X(ceil(win/2):end-floor(win/2),ceil(win/2):end-floor(win/2),:);
  [m,n,p] = size(X);
  X = reshape(X,m*n,p)';
  Index = ceil(find(X>MaxPurity)/p);
  X(:,Index) = 1/p*ones(p,1)*ones(1,length(Index));
  outliers=10;
  violation_extremes = [-1 0];
  spread = violation_extremes(2)-violation_extremes(1);
for i=1:outliers
    index= randperm(p);
    X(index(1),i) = violation_extremes(1) + spread*rand(1);
    aux = rand(p-1,1);
    aux = aux./repmat(sum(aux),p-1,1)-X(index(1),i)/(p-1);
    X(index(2:p),i) = aux;
end
   HIM = reshape((M*X)',m,n,bands);    
   Y = reshape(HIM,[m*n,bands,1]);
   SNR=30;
 if SNR <100
% add noise
   variance = sum(Y(:).^2)/10^(SNR/10)/m/m/bands;
   n = sqrt(variance)*randn([bands m*n]);
   Y = Y' + n;
 else
   Y=Y';
end
% display
[y1,Up,~,~] = dataProj(Y,p,'proj_type','affine');
y1=Up'*y1;
Mtrue = Up'*M;
I = 1;
J = 2;
K = 3;
E_I = eye(p);
v1 = E_I(:,I);
v2 = E_I(:,J);
v3 = E_I(:,K);
y1 = [v1 v2 v3]'*y1;
m_true = [v1 v2 v3]'*Mtrue;
mm = m_true;
figure
scatter3(y1(1,:),y1(2,:),y1(3,:),20);
hold on;
plot3(m_true(1,[1:end 1]), m_true(2,[1:end 1]),m_true(3,[1:end 1]),'black-','MarkerSize',30);
plot3(m_true(1,[1:end 2]),m_true(2,[1:end 2]),m_true(3,[1:end 2]), 'black-','MarkerSize',30);
plot3(m_true(1,[1:1 3]),m_true(2,[1:1 3]),m_true(3,[1:1 3]), 'black-','MarkerSize',30);
plot3(mm(1,:),mm(2,:),mm(3,:), 'black.','MarkerSize',10);
scatter3(y1(1,1:outliers),y1(2,1:outliers),y1(3,1:outliers),80,'rx');
% title('Linear Model','fontsize',18);
end