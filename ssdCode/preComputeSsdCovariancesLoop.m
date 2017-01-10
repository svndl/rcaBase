function [sum_S,sum_N,nPointsIn_S,nPointsIn_N]=preComputeSsdCovariancesLoop(data_signal,data_noise1,data_noise2,condRange,subjRange)
% [SUMXX,SUMYY,SUMXY,NPOINTSINXX,NPOINTSINYY,NPOINTSINXY]=PRECOMPUTERCACOVARIANCESLOOP(DATA)
%
% take in epoched data volume(s) to produce signal and noise covariance matrices for input into ssdTrain()
%
% data can either be a 3-D array (samples x channels x trials), or a cell
% array (conditions x subjects) of 3-D data arrays (as above)
%
% throughout data, NaNs indicate missing data values (i.e., artifacts)
%
% sum_S, sum_N: unnormalized auto-covariance of the signal and the noise, repectively
% nPointsIn_S, nPointsInN: number of non-missing data points used to compute covariances
%
% covariances are left unnormalized to allow for subsequent averaging
% across subjects or conditions, as appropriate
%
% file is an adaptation of preComputeRcaCovariancesLoop.m by Jacek P. Dmochowski
% SB, 2016

assert(sum(size(data_signal) == size(data_noise1) & size(data_signal) == size(data_noise2))==length(size(data_signal)));
if ~iscell(data_signal),
    nCond=1; nSubjects=1; subjRange=1; condRange=1; data_signal={data_signal}; data_noise1={data_noise1}; data_noise2={data_noise2};
else
    if nargin<5, subjRange=1:size(data_signal,2); end
    if nargin<3, condRange=1:size(data_signal,1); end
    data_signal=data_signal(condRange,subjRange);data_noise1=data_noise1(condRange,subjRange);data_noise2=data_noise2(condRange,subjRange);
    nCond=size(data_signal,1);
    nSubjects=size(data_signal,2);
end

for subj_idx = 1:nSubjects
    for cond_idx = 1:nCond
        assert(sum(size(data_signal{cond_idx,subj_idx}) == size(data_noise1{cond_idx,subj_idx}) & size(data_signal{cond_idx,subj_idx}) == size(data_noise2{cond_idx,subj_idx}))==length(size(data_signal{cond_idx,subj_idx})));
    end
end

fprintf('Selected %d subjects and %d condition(s) for training... \n',nSubjects,nCond);

[~,nElectrodes,~] = size(data_signal{1,1});  % assume uniform (?)
sum_S             = zeros(nCond,nSubjects,nElectrodes,nElectrodes); sum_N=zeros(nCond,nSubjects,nElectrodes,nElectrodes);
nPointsIn_S       = zeros(nCond,nSubjects,nElectrodes,nElectrodes); nPointsIn_N=zeros(nCond,nSubjects,nElectrodes,nElectrodes);

%if verLessThan('matlab','8.2') % 8.2 = version number for R2013b
%    matlabpool;
%else
%    parpool; % works for all versions of matlab from R2013b forward (matlabpool was removed in R2015a)
%end

for cond=1:nCond
    for subj=1:nSubjects
        fprintf('Computing covariances for subject %d/%d and condition %d/%d... \n',subjRange(subj),nSubjects,condRange(cond),nCond);
        thisVolume_signal=permute(data_signal{cond,subj},[3,2,1]); % trials, channels, frequencies * bins
        thisVolume_noise1=permute(data_noise1{cond,subj},[3,2,1]);
        thisVolume_noise2=permute(data_noise2{cond,subj},[3,2,1]);
        [nTrials,~,nSamples]=size(thisVolume_signal);
        if nSamples<nElectrodes, warning('Number of samples is less than the number of electrodes'); end
        
        C_s = zeros(nElectrodes); C_n = zeros(nElectrodes);
        P_s = zeros(nElectrodes); P_n = zeros(nElectrodes);
        
        if nSamples == 1 
            % check if more than one (harmonic) frequency is contained in X_s.
            % Due to use of transpose() and squeeze() functions in covariance
            % matrix computation, the case in which there is only 1 frequency has
            % to be treated differently than cases in where there are multiple
            % frequencies
            

            for e = 1:nTrials % e stands for epoch
                error('there are still errors in the branch!') % SB TODO
                % Parseval's theorem leads to the expression
                % 2*real(sum_freq(a * conj(b))). Due to the SSVEP, we only
                % expect harmonics of the stimulation frequency to carry
                % any energy, so the sum is only accross harmonics.
                % Furthermore, we use the relation 2*real(z) = z + conj(z). 
                
                M_S             =           [transpose(thisVolume_signal(e,:)), transpose(conj(thisVolume_signal(e,:)))];
                M_N1            =           [transpose(thisVolume_noise1(e,:)), transpose(conj(thisVolume_noise1(e,:)))];
                
                
                M_S_transposed  = transpose([transpose(conj(thisVolume_signal(e,:))), transpose(thisVolume_signal(e,:))]);
                % M_N1, M_N1_t, M_N2, M_N2_t are analogous to M_S and M_S_t
                
                M_N1_transposed = transpose([transpose(conj(thisVolume_noise1(e,:))), transpose(thisVolume_noise1(e,:))]);
                M_N2            =           [transpose(thisVolume_noise2(e,:)), transpose(conj(thisVolume_noise2(e,:)))];
                M_N2_transposed = transpose([transpose(conj(thisVolume_noise2(e,:))), transpose(thisVolume_noise2(e,:))]);
                
                % Construct covariance matrix by summing over covariances
                % from each trial
                C_s = C_s+ M_S*M_S_transposed;
                C_n = C_n+ M_N1*M_N1_transposed + M_N2*M_N2_transposed;
                
                P_s = P_s + double(~isnan(M_S))*double(~isnan(M_S_transposed));
                P_n = P_n + double(~isnan(M_N1)*double(~isnan(M_N1_transposed))) +double(~isnan(M_N2))*double(~isnan(M_N2_transposed));
            end
            
        else
            % same as above but with multiple frequencies, i.e., multiple
            % harmonics
            
            % vectorized computation of trial averaged covariance matrix
            %C_s = zeros(nElectrodes); C_n = zeros(nElectrodes);
            %P_s = zeros(nElectrodes); P_n = zeros(nElectrodes);
            for e = 1:nTrials % loop over epochs -> trial-averaged Covariance matrix
%                 
%                 % thisVolume_signal: trials, channels, frequencies*bins
%                 M_S             = [squeeze(thisVolume_signal(e,:,:)), conj(squeeze(thisVolume_signal(e,:,:)))];
%                 M_N1            = [squeeze(thisVolume_noise1(e,:,:)), conj(squeeze(thisVolume_noise1(e,:,:)))];
%                 M_N2            = [squeeze(thisVolume_noise2(e,:,:)), conj(squeeze(thisVolume_noise2(e,:,:)))];
%                             
%                 M_S_transposed  = transpose([conj(squeeze(thisVolume_signal(e,:,:))), squeeze(thisVolume_signal(e,:,:))]);
%                 M_N1_transposed = transpose([conj(squeeze(thisVolume_noise1(e,:,:))), squeeze(thisVolume_noise1(e,:,:))]);
%                 M_N2_transposed = transpose([conj(squeeze(thisVolume_noise2(e,:,:))), squeeze(thisVolume_noise2(e,:,:))]);
%                 
%                 P_s = P_s + double(~isnan(M_S))*double(~isnan(M_S_transposed));
%                 P_n = P_n + double(~isnan(M_N1+M_N2))*double(~isnan(M_N1_transposed+M_N2_transposed)) ;
%                 
%                 M_S(isnan(M_S)) = 0; M_S_transposed(isnan(M_S_transposed)) = 0;
%                 M_N1(isnan(M_N1)) = 0; M_N1_transposed(isnan(M_N1_transposed)) = 0;
%                 M_N2(isnan(M_N2)) = 0; M_N2_transposed(isnan(M_N2_transposed)) = 0;
%                 
%                 %thisTrialCovMat_S = [squeeze(thisVolume_signal(e,:,:)), conj(squeeze(thisVolume_signal(e,:,:)))] * transpose([conj(squeeze(thisVolume_signal(e,:,:))), squeeze(thisVolume_signal(e,:,:))]);
%                 %thisTrialCovMat_N = [squeeze(thisVolume_noise1(e,:,:)), conj(squeeze(thisVolume_noise1(e,:,:)))] * transpose([conj(squeeze(thisVolume_noise1(e,:,:))), squeeze(thisVolume_noise1(e,:,:))])+...
%                 %                    [squeeze(thisVolume_noise2(e,:,:)), conj(squeeze(thisVolume_noise2(e,:,:)))] * transpose([conj(squeeze(thisVolume_noise2(e,:,:))), squeeze(thisVolume_noise2(e,:,:))]);
%                 
%                 C_s = C_s + M_S*M_S_transposed;
%                 C_n = C_n + M_N1*M_N1_transposed + M_N2 * M_N2_transposed;
                temp_s  = squeeze(thisVolume_signal(e,:,:)) ;
                temp_n1 = squeeze(thisVolume_noise1(e,:,:)) ;
                temp_n2 = squeeze(thisVolume_noise2(e,:,:)) ;
                
                P_s = P_s + double(~isnan(temp_s))*double(~isnan(temp_s))' ; 
                P_n = P_n + 0.5*double(~isnan(temp_n1))*double(~isnan(temp_n1))'; 
                P_n = P_n + 0.5*double(~isnan(temp_n2))*double(~isnan(temp_n2))' ;
                
                temp_s(isnan(temp_s)) = 0 ;
                temp_n1(isnan(temp_n1)) = 0 ;
                temp_n2(isnan(temp_n2)) = 0 ;
                
                C_s  = C_s + [temp_s,conj(temp_s)]*[temp_s,conj(temp_s)]';
                C_n  = C_n + [temp_n1,conj(temp_n1)]*[temp_n1,conj(temp_n1)]';
                C_n  = C_n + [temp_n2,conj(temp_n2)]*[temp_n2,conj(temp_n2)]';
            end
        end
        sum_S(cond,subj,:,:) = C_s; sum_N(cond,subj,:,:) = C_n;
        nPointsIn_S(cond,subj,:,:) = P_s; nPointsIn_N(cond,subj,:,:) = P_n;
    end
end

sum_S=squeeze(sum_S); sum_N=squeeze(sum_N);
if max(abs(imag(sum_S(:)))) > 1e-10
    warning('imaginary values in S!')
end

if max(abs(imag(sum_N(:)))) > 1e-10
    warning('imaginary values in S_n!')
end
sum_S=real(sum_S);
sum_N=real(sum_N);
    
nPointsIn_S=squeeze(nPointsIn_S); nPointsIn_N=squeeze(nPointsIn_N);

%if verLessThan('matlab','8.2')
%    matlabpool close;
%else
%    delete(gcp); %closes the parpool
%end
