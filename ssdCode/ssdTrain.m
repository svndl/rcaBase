% define the ssd function for our application based on signal and noise
% definitions in the fourier domain
function [Wsub, A, W, dGenSort, K] = ssdTrain(C_s, C_n, nComp, K)

%     X_s: matrix containing the complex fourier coefficients at the frequency of interest extracted  epoch-wise
%     X_n2, X_n2: matrix containing the complex fourier coefficients at the neighbor bins of the  frequency of interest at different epochs
%     note that in this implementation, coefficients are only given for positive frequencies. parseval's theorem is met by
%     considering coefficients at negative fequencies as the conjugates of coefficients at positive frequencies

%   Soren:
%   X_s, X_n1, X_n2: dimensions are [epochs x channels x frequencies]
%   (i.e. [trials x channels x frequencies]).

assert(~any(~(size(C_s) == size(C_n))), 'X_s, X_n1, X_n2 should have the same dimensions.');

if (mean(abs(imag(C_s(:))))>10^-10 ) | (mean(abs(imag(C_n(:))))>10^-10 )
    warning('something went wrong!')
end

if sum(isnan(C_s(:)))~=0 || sum(isnan(C_n(:)))~=0
    warning('Covariance matrices contain NaNs: setting to zero...');
    C_s(isnan(C_s))=0; C_n(isnan(C_n))=0;
end

C_s = real(C_s); % When tested, max(imag(C_s(:))) was in the order of e-17. It is discarded as a numerical artifact.
C_n = real(C_n);

% regurlarized based on covariance of signal
[Vpool,Dpool]=eig(C_s);
[dPool,sortIndx]=sort(diag(Dpool),'ascend');
Vpool=Vpool(:,sortIndx);  % in case eigenvectors/eigenvalues not sorted
dPool=dPool(end-K:end);
Rw = C_n \ Vpool(:,end-K:end)*diag(dPool)*Vpool(:,end-K:end)';


[Vgen,Dgen]=eig(Rw);  % compute generalized eigenvalues/eigenvectors
dGen=diag(Dgen);
[dGenSort,b]=sort(abs(dGen));  
W=Vgen(:,b(end:-1:1)); % in case not sorted

if max(abs(imag(W))) > 10e-10
    warning('W has imaginary part')
end

W=real(W);  % ignore small imaginary component

Wsub=W(:,1:nComp);  % return only selected number of components
A=C_s*Wsub / (Wsub'*C_s*Wsub);  % return forward models (see Parra et al., 2005, Neuroimage)

end
