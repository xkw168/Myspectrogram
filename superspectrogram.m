%______________________________________________________________________________________________________________________
% 
% @file         superspectrogram.m
% @date         27/11/2017
% @author       Kewei Xia. Xinhai Pan
% @affiliation  2017 Digital Signal Processing Course Design, 15 EE
% @brief        myspectrogram function
%______________________________________________________________________________________________________________________
%
% @inputs       speech  - time domain speech signal vector
%               fs      - sampling frequency (Hz), f.e. 8000
%               Fr      - Display Frequency range vector (Hz), f.e. [0 2000] 
%               T       - vector of frame width, Tw, and frame shift, Ts, in milliseconds, i.e. [Tw, Ts]
%               w       - analysis window choice, f.e. @hamming, @hanning
%               nfft    - fft analysis length, f.e. 1024
%               Slim    - vector of spectrogram limits (dB), i.e. [Smin Smax]
%               alpha   - fir pre-emphasis filter coefficients, f.e. [1 -0.97]
%               cmap    - color map ('default', 'gray', 'bone', 'copper', 'hot', 'jet')
%               cbar    - color bar (boolean)
%
% @output       handle  - plot handle
%______________________________________________________________________________________________________________________
%
% @usage        [handle] = superspectrogram(s, fs, Fr, T, w, nfft, Slim, alpha, cmap, cbar);
% @examples     [handle] = superspectrogram(speech, 8000,[0 2000], [18 1], @hamming, 1024, [-45 -2], false, 'default', false);
%               [handle] = superspectrogram(speech, 8000,[0 4000], [18 1], @hanning, 512, [-50 -2], [1 -0.97], 'default', true);
%______________________________________________________________________________________________________________________
function [handle] = superspectrogram(s, fs, Fr, T, w, nfft, Slim, alpha, cmap, cbar)
    %__________________________________________________________________________________________________________________
    % verify input and set defaults
    switch nargin
            case 1, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[18,1]; nfft=1024; fs=8000; Fr=[0,0];
            case 2, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[18,1]; nfft=1024; Fr=[0,0];
            case 3, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[18,1]; nfft=1024;
            case 4, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; nfft=1024;
            case 5, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; nfft=1024;
            case 6, cbar=true; cmap='default'; alpha=false; Slim=[-59,-1];
            case 7, cbar=true; cmap='default'; alpha=false;
            case 8, cbar=true; cmap='default';
            case 9, cbar=true;
            case 10
        otherwise
            error('Invalid number of input arguments.');
    end
    
    %__________________________________________________________________________________________________________________
    % declare variables
    if(ischar(s))
        [s, fs] = audioread(s);                                  % read audio data from file
    end
	
	fmin = Fr(1);
    fmax = Fr(2);
	if(fmin == fmax)
        fmax = fs/2;
    end
	
    Tw = T(1);                                                        % frame width (ms)
    Ts = T(2);                                                         % frame shift (ms)
    Nw = round(fs*Tw*0.001);                               % frame width (samples) 
    Ns = round(fs*Ts*0.001);                                 % frame shift (samples) 
    N  = length(s);                                                 % length of speech signal (samples)
    Smin = Slim(1);                                                % lower normalized dynamic range limit 
    Smax = Slim(2);                                               % upper normalized dynamic range limit
    if(ischar(w))
        w = str2func(w);                                          % obtain window function handle from string input
    end 

    %__________________________________________________________________________________________________________________
    % preprocessing input speech signal
    if(islogical(alpha) && alpha)
        s = filter([1 -0.95],1,s);                                 % apply a typical preemphasis filter(FIR)
    elseif(~islogical(alpha)) 
        s = filter(alpha,1,s);                                      % apply custom preemphasis filter, by default
    end       

    %__________________________________________________________________________________________________________________
    % get spectrogram data, using function: toframes, which will be declared later 
    [S,F,T] = toframes(s, [fmin,fmax], w, T, fs, nfft);         % Framing function

    %__________________________________________________________________________________________________________________
    % set dynamic spectrogram range 
    S = abs(S);                                                       % compute magnitude spectrum 
    S = S/max(max(S));                                          % normalize magntide spectrum
    S = 20*log10(S);                                               % compute power spectrum in dB
    
    F = F(round(fmin/fs*nfft)+1 : round(fmax/fs*nfft));
    
    %__________________________________________________________________________________________________________________
    % plot using function: imagesc 
    handle = imagesc(T, F, S, [Smin Smax]);
    axis('xy');
    axis([0 N/fs  round(fmin) round(fmax)]);
	%axis([0 N/fs  0 fs/2]);
    xlabel('time (s)', 'FontSize', 8, 'FontWeight', 'n');
    ylabel('frequency (Hz)', 'FontSize', 8, 'FontWeight', 'n');
    set(gca,'YDir','normal', 'FontSize', 6);

    if(cbar)
        colorbar('FontSize',6); 
    end

    %__________________________________________________________________________________________________________________
    % define default colormap: reverse gray colormap for better display
    switch(lower(cmap))
    case {'default'}
        colormap('gray');
        map=colormap;
        colormap(1-map);    
    otherwise
    	colormap(cmap);
    end

%______________________________________________________________________________________________________________________
% 
% @author       Kewei Xia, Xinhai Pan
% @date         11, 2017
% @brief        framing function to get spectrogram data
%______________________________________________________________________________________________________________________
function [S,F,T] = toframes(s, Fr, w, T, fs, nfft)

    %__________________________________________________________________________________________________________________
    % verify inputs
    switch nargin
        case 1, nfft=1024; fs=8000; T=[32 32/8]; w=@hamming; Fr = [0,0];
        case 2, nfft=1024; fs=8000; T=[32 32/8]; w=@hamming;
        case 3, nfft=1024; fs=8000; T=[32 32/8];
        case 4, nfft=1024; fs=8000;
        case 5, nfft=1024;
        case 6
    otherwise
    	error('Invalid number of input arguments.');
    end

    %__________________________________________________________________________________________________________________
    % define variables
    if(ischar(s)) 
        [s, fs] = audioread(s);
    end

    s = s(:).';
    smax = max(abs(s));
    s = s/smax;                                                      % normalize

	fmin = Fr(1);
    fmax = Fr(2);
	if(fmin == fmax)
	    fmax = fs / 2;
	end
	
    Tw = T(1);                                                        % frame length [ms]
    Ts = T(2);                                                         % frame frameshift [ms]
    N  = round(fs*Tw*0.001);                                % frame length [elements]
    Z  = round(fs*Ts*0.001);                                  % frame shift [elements]
	D  = mod(length(s-N), Z);                               % add N-D zeros to the end, padding
	s  = [s zeros(1,N-D)];                                      % assure the number of frame is integer
    ss = length(s);                                                 % length of the signal for processing
    M  = ((ss-N)/Z)+1;                                          % number of overlapping segments
    
    %__________________________________________________________________________________________________________________
    % identify window type function
    if(strcmp(w,'@hamming'))
        w=str2func(w);
        wa = w(N).';
    elseif(strcmp(w,'@hanning'))
        p=str2func(w);
        wa = p(N).';
    elseif(strcmp(w,'@blackman'))
        q=str2func(w);
        wa = q(N).';
    end
    wa = w(N).';                                                     % default


    %__________________________________________________________________________________________________________________
    % split the signal into frames and windowing
    indf = Z*(0:(M-1)).';                                         % indexes for frames
    inds = (1:N);                                                    % indexes for samples
    refs = indf(:,ones(1,N)) + inds(ones(M,1),:);    % sample indexes for each frame
    
    segments_s = s(refs);                                      % split into overlapped frames (using indexing)
    segments_sm = segments_s .* wa(ones(M,1),:);    % apply magnitude spectrum analysis window 


    %__________________________________________________________________________________________________________________
    % apply fft to segments_sm
    F = [0:nfft-1]/(nfft-1)*fs;                                  %0-fs the sampling rate of the signal
    T = [0:M-1]/(M-1)*ss/fs;                                  %0-t  the length of the audio
    S = fft(segments_sm, nfft, 2);                          % short-time Fourier analysis for each row
	

	S = S(:,[round(fmin/fs*nfft)+1 : round(fmax/fs*nfft)]); %extract the desire frequency range
    S = abs(S).^2/N;                                             % periodogram PSD estimates
    S = sqrt(S);                                                      % magnitude spectrum (for consistency)
    S = S.';
    F = F.';
    T = T;

%______________________________________________________________________________________________________________________
% 
% @author       Kewei Xia, Xinhai Pan
% @date         Nov, 2017
% @brief        Hamming window function
%______________________________________________________________________________________________________________________
function [w] = hamming(N)

    w = 0.54-0.46*cos(2*pi*[0:(N-1)].'/(N-1));

%______________________________________________________________________________________________________________________
% EOF


%______________________________________________________________________________________________________________________
% 
% @author       Kewei Xia, Xinhai Pan
% @date         Nov, 2017
% @brief        Hanning window function
%______________________________________________________________________________________________________________________
function [p] = hanning(N)

   p=0.5*(1+cos(2*pi*[0:(N-1)].'/(N-1)));
   
%______________________________________________________________________________________________________________________
% EOF


%______________________________________________________________________________________________________________________
% 
% @author       Kewei Xia, Xinhai Pan
% @date         Nov, 2017
% @brief        Hanning window function
%______________________________________________________________________________________________________________________
function [q] = blackman(N)

    q = 0.42-0.5*cos(2*pi*[0:(N-1)].'/(N-1))+0.08*cos(4*pi*[0:(N-1)].'/(N-1));
    
%______________________________________________________________________________________________________________________
% EOF
