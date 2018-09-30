%some signal for test

%chirp signal
fs=8000;
ts=1/fs;
time=0:ts:2;
freqs=[500 590];
s=zeros(length(freqs),length(time));
for i=1:length(freqs)
    xs(i,:)=cos(2 * pi * freqs(i) * time + 200 * cos(2 * pi * time));
end
x=sum(xs);
x=x./max(abs(x));

% single-frequency signal
Fs = 8000; % 采样频率
T = 4;     % 时间长度
n = Fs * T;  % 采样点数
f = 500;   % 声音频率Hz
y = sin(2 * pi * f * T * linspace(0, 1, n + 1));
sound(y, Fs);








