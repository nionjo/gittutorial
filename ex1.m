close all
clear all

load binarydigits.txt -ascii;
Y=binarydigits;
[N D]=size(Y);
maxiter=100;

Lfinal=zeros(3,1)
Pi=zeros
for k=2:4
    [Pi,P,l, R, lprev, iter, L]=EM(k,Y,maxiter);
    figure; plot(L);
    title(sprintf('K = %i', k))
    Lfinal(k-1)=log(abs(L(end)))
    P_{k-1}=P;
    Pi_{k-1}=Pi;
    R_{k-1}=R
    figure
    colormap gray;
    for j=1:k
        subplot(1,k,j);
        imagesc(reshape(P(j,:)',8,8)'); % plot mixing components
        axis off;
    end
end


