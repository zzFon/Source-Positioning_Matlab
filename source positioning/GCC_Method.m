function G = GCC_Method(m,s1,s2,wnd,inc)
    % ���廥���ʱ�ӹ���
    % ��2��mic���յ����źŽ��й��廥��أ��Ӷ�����2��mic���ʱ�ӹ���
    
    % ���廥��صļ�Ȩ�������ͣ�'standard' 'roth' 'scot' 'phat' 'ml'
    % mic1���յ��ź�s1
    % mic2���յ��ź�s2
    % ���ڳ���wnd
    % ֡��inc
    % 2��mic�����ʱ��G
    
    N = wnd;
    wnd = hamming(N); % (4-3)hamming��
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
