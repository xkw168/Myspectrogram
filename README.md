# myspectrogram
myspetrogram（）的代码，功能十分齐全，运行环境matlab2017a

	
	@file         read me.txt
	@date         27/11/2017
	@author       Kewei Xia. Xinhai Pan
	@affiliation  2017 Digital Signal Processing Course Design, 15 EE
	@brief        instructions
	
##编译环境
MATLAB R2016a		注：版本过低可能会导致运行问题，主要为函数变更。例如：wavread——audioread
##文件注释
1. signal.m: 用于生成测试信号。x：chirp signal   y：single freq signal
2. superspectrogram.m: 函数文件，不直接运行；

##测试流程
利用MATLAB软件打开signal.m 文件，打开编译器界面点击运行即可得出测试结果。
要查看函数代码，则利用MATLAB软件打开superspectrogram.m文件。

---
###superspectrogram.m	
>superspectrogram(s, fs, Fr, T, w, nfft, Slim, alpha, cmap, cbar);
%用户交互接口，调用方式见注释及signal.m
toframes(s, Fr, w, T, fs, nfft);
%信号截断/处理函数，默认窗长为32，默认帧移为4的重叠窗
hamming(N);
%hamming窗函数，返回N点hamming窗
hanning(N);
%hanning窗函数，返回N点hanning窗
blackman(N);
%blackman窗函数，返回N点blackman窗
%注：对于带限信号，一般采用以上三角窗进行处理，本实验首先装载了
%三种窗函数，其他窗型有待加载，例如kaiser窗

###signal.m
>x: chirp signal 其stft为正弦状；	y: single-frequency signal 单频信号
figure(1): superspectrogram display: in default conditions
%默认参数情况下调用显示x,y的语谱图
figure(2): superspectrogram display: adjust spectrogram limits
%改变显示频率范围/显示衰减范围，提高显示效果
figure(3): superspectrogram display: adjust window type
%改变窗型，对比显示效果


