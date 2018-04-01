function [H] = SNSAabun(y_augment, A_augment, P, MaxIteration, Threshold);
W = A_augment;
inp = y_augment;
H = P;
alpha = 1; beta = 0.1; sigma = 0.01;
for iteration = 1: MaxIteration
    tW = W'*W*H - W'*inp;
    gradient = norm(tW(tW < 0 | H >0));
    if gradient <  Threshold,
        break
    end
% gradient descent
    for initeration=1:50,
        Hn = max(H - alpha*tW, 0); 
        dN = Hn-H;
        tWW=sum(sum(tW.*dN)); 
        dd = sum(sum((W'*W*dN).*dN));
        suff_decr = 0.99*tWW + 0.5*dd < 0;
        if initeration == 1
            decr_alpha = ~suff_decr;
            Hh = H;
        end
        if decr_alpha 
            if suff_decr,
                H = Hn; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr | Hh == Hn,
                H = Hh; break;
            else
                alpha = alpha/beta; Hh = Hn;
            end
        end
    end
end
end