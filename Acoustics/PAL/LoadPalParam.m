function [transducer, ultra, audio] = LoadParam(varargin)

    p = inputParser;
    addParameter(p, 'set_name', 'ZhongNew');
    addParameter(p, 'AudioFreq', 1e3);
    parse(p, varargin{:});
    ip = p.Results;

    audio = SoundWave(ip.AudioFreq);

    switch ip.set_name
        % the parameters from the paper: M. Červenka and M. Bednařík, ...
        %   “A versatile computational approach for the numerical ...
        %   modelling of parametric acoustic array,” J. Acoust. ...
        %   Soc. Am. 146(4), 2163-2169 (2019).
        case 'Cervenka'
            ultra.low = SoundWave(39e3);
            ultra.high = SoundWave(ultra.low.freq + audio.freq);

            % 超声振速
            ultra.vel0_low = 1;
            ultra.vel0_high = 1;

            delta = 3.764e-5;
            alpha = ultra.low.angularFreq^2*delta /2/343^3;
            ultra.low.wavnumber = ultra.low.wavnumber + 1i*alpha;
            ultra.high.wavnumber = ultra.high.wavnumber + 1i*alpha;
                    
            % 换能器参数
            transducer.radius = 0.02;
        case 'Zhong'
            ultra.low = SoundWave(64e3);
            ultra.high = SoundWave(ultra.low.freq + audio.freq);
            ultra.vel0_low = 0.12;
            ultra.vel0_high = 0.12;
            
            % 换能器参数
            transducer.radius = 0.1;

        case 'ZhongNew'
            ultra.ave = SoundWave(64e3);
            ultra.low = SoundWave(ultra.ave.freq - audio.freq/2);
            ultra.high = SoundWave(ultra.ave.freq + audio.freq/2);
            ultra.vel0_low = 0.12;
            ultra.vel0_high = 0.12;
    end
end
