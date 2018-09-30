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
Fs = 8000; % ����Ƶ��
T = 4;     % ʱ�䳤��
n = Fs * T;  % ��������
f = 500;   % ����Ƶ��Hz
y = sin(2 * pi * f * T * linspace(0, 1, n + 1));
sound(y, Fs);








