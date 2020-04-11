function G = GCC_Method(m,s1,s2,wnd,inc)
    % 广义互相关时延估计
    % 用2个mic接收到的信号进行广义互相关，从而进行2个mic间的时延估计
    
    % 广义互相关的加权函数类型：'standard' 'roth' 'scot' 'phat' 'ml'
    % mic1接收到信号s1
    % mic2接收到信号s2
    % 窗口长度wnd
    % 帧移inc
    % 2个mic的相对时延G
    
    N = wnd;
    wnd = hamming(N); % (4-3)hamming窗
    x = enframe(s1,wnd,inc);
    y = enframe(s2,wnd,inc);
    n_frame = size(x,1);
    
    switch lower(m)
        case 'standard'
            for i = 1:n_frame
                x = s1(i:i+N)
                y = s2(i:i+N);
                X = fft(x,2*N-1);
                Y = fft(y,2*N-1);
                Sxy = X.*conj(Y);
                gain = 1;
                Cxy = fftshift(ifft(Sxy.*gain));
                [Gvalue(i),G(i)] = max(Cxy);
            end
        otherwise error('Method not defined');
end
