function [U,P] = NFINDR_Sun(HIM,p);
% NFINDR algorithm for endmember extraction.
% ------------------------------------------------------------------------------
% Input:   HIM : hyperspectral image cube [nrows x ncols x nchannels]
%          p   : desired number of endmembers to be extracted
%
% Output:  U   : set of extracted endmembers [nchannels x p] 
%          P   : spatial coordinates of the extracted endmembers (positions
%                rows x cols)
%
% Note 1 -     in this NFINDER version, the PCA transformation is used for 
%              the initial (p-1) feature reduction, and the extracted features
%              are used to calculate the volume.
%
% Note 2 -     once the first NFINDR solution is found (original idea),
%              we run again the optimization process to try a further improvement
%              in the final volume.
%
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 

KindFeatRed = 'pca';

disp(' === Start NFINDR run ===')

if strcmp(KindFeatRed,'pca') % USES THE PCA FOR  FEATURE REDUCTION

    disp('. this NFINDR version uses the PCA for the iniciatl (p-1) feature reduction')

    % Obtener tiempo CPU actual
    t1=cputime;

    [nr,nc,nb] = size(HIM);
    aux = reshape(HIM,[nr*nc,nb,1]);

    [pc, zscores, pcvars] = princomp(aux);
    cc = cumsum(pcvars./sum(pcvars) * 100);
    
    disp(sprintf('. Variances:'))
    disp(sprintf(' %4.2f ',pcvars'))
    disp(sprintf('. Cumulate variance of the first [%d] components = %4.2f /100',p-1,cc(p-1)))

    MNFHIM = reshape(zscores,[nr,nc,nb]);

else
    error(' please provide a valid feature reduction method...')
end

% Obtener tiempo CPU actual
t2=cputime;

[ns,nl,nb] = size(HIM);
P = zeros(2,p);

MatrixInicial = [];
rand('seed',0);
for k = 1:p
    row = floor(rand*ns+1);
    line = floor(rand*nl+1);
    disp(['Random pixel in position:',num2str(k),':',num2str(row),',',num2str(line)]);
    MatrixInicial = [ MatrixInicial squeeze(MNFHIM(row,line,1:p-1)) ];
    P(1,k) = row; P(2,k) = line;
end

MatrixTest = zeros(p,p);
MatrixTest(1,:) = 1;
MatrixTest(2:p,:) = MatrixInicial;

volumeactual = abs(det(MatrixTest)); % instead of: volumeactual = abs(det(MatrixTest))/(factorial(p-1));

maxit = 3*p; % heuristic  in this version, allowing an additional optimization! 
%disp(sprintf('. Maximun # of iterations allowed [%d]x[%d]=[%d] (heuristic inplemented in this version)',3,p,maxit))

it = 1;
v1 = -1;
v2 = volumeactual;

t4=cputime;
while and(it<=maxit,v2>v1)

    disp(sprintf('. Start iteration # [%d], with abs(det(E)) = %4.8g',it,v2))
    P
    disp(sprintf('  Functional abs(det(E)): @ v1 = %5.5g  @ v2 = %5.5g , Ratio: v2/v1 = %5.5g ',v1,v2,v2/v1))

    for k = 1:p

        disp(sprintf('. Loop @ endmember # %d',k))

        for i = 1:ns
            for j = 1:nl
                actual = MatrixTest(2:p,k);
                MatrixTest(2:p,k) = squeeze(MNFHIM(i,j,1:p-1));
                volume = abs(det(MatrixTest));  % instead of: volume = abs(det(MatrixTest))/(factorial(p-1));
                if volume > volumeactual
                    disp(sprintf('---> update with pixel @ (%5d,%5d) | abs(det(E)) = %5.5g ',i,j,volume))
                    volumeactual = volume;
                    P(1,k) = i; 
                    P(2,k) = j;
                else
                    MatrixTest(2:p,k) = actual;
                end
            end
        end
    end
    
    disp(sprintf('. End of iteration # [%d]: abs(det(E)) = %5.5g, the pixels are:',it,volumeactual))
    P
    
    it = it+1;

    v1 = v2;
    v2 = volumeactual;

end
t5=cputime;

valfunct = volumeactual;
nit = it-1;

if nit<maxit
    disp(sprintf(' End, convergence @ iteration # [%d]. Final abs(det(E)) = %5.5g',nit,volumeactual))
else
    disp(sprintf(' End, NO convergence until iteration # %d. The abs(det(E)) = %5.5g',nit,volumeactual))
end

disp(sprintf('. The NFINDER found solution is:'))
P = P'

U = [];
for i = 1:p
    U = [ U squeeze(HIM(P(i,1),P(i,2),:)) ];
end

% Visualizar HIMn
% figure
% imagesc(mean(HIM,3)); colormap(gray);
% set(gca,'DefaultTextColor','black','xtick',[],'ytick',[],'dataaspectratio',[1 1 1])

% for i=1:size(P,1)
%     % Visualizar pixel seleccionado en pantalla
% drawnow;
% %    text(P(i,2),P(i,1),'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','green');
% text(P(i,2),P(i,1),'O','Color','red');
%  
%  end

% Obtener tiempo CPU actual
t3=cputime;

% Mostrar tiempo total en ejecucion del algoritmo
disp(sprintf('. CPU Processing time for feature reduction..... %6.3f [s]  ',(t2-t1)));
disp(sprintf('. CPU Processing time for iterative loop   ..... %6.3f [s]  ',(t5-t4)));
disp(sprintf('. Total CPU processing time .................... %6.3f [s]  ',(t3-t1)));
disp(' === Eng NFINDR ===');


