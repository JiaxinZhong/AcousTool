% =========================================================================
%	Provide the liset of octave frequencies and 1/3 octave frequencies
% -------------------------------------------------------------------------
% Last modified --- 23-Jan-2019
% =========================================================================
function center_freq = Octave(freq_low, freq_up)

    freq_list = [25; 31.5; 40; 50; 63; 80; 100; 125; 160; 200; 250; 315; 400; 500; 630; 800; 1000; 1250; 1600; 2000; 2500; 3150; 4000; 5000; 6300; 8000; 10000; 12500; 16000; 20000];
    
    center_freq = freq_list(freq_list >= freq_low & freq_list <= freq_up);

end

% classdef Octave < handle
% 
% 	properties
% 		CENTER_FREQUENCY
% 		CENTER_FREQUENCY_3
%         center_freq
%         center_freq3
% 	end
% 
% 	methods 
% 		function obj = Octave(freqStart, freqEnd)
% 			% Center frequencies for 1/3 octave
% 			obj.CENTER_FREQUENCY_3 = [25; 31.5; 40; 50; 63; 80; 100; 125; 160; 200; 250; 315; 400; 500; 630; 800; 1000; 1250; 1600; 2000; 2500; 3150; 4000; 5000; 6300; 8000; 10000; 12500; 16000; 20000];
% 			obj.center_freq3 = ...
%                 obj.CENTER_FREQUENCY_3(obj.CENTER_FREQUENCY_3 >= freqStart & ...
%                 obj.CENTER_FREQUENCY_3 <= freqEnd);
%             
% 			% Center frequencies for octave
% 			obj.CENTER_FREQUENCY = obj.CENTER_FREQUENCY_3(2 : 2 : end);
%             obj.center_freq = ...
%                 obj.CENTER_FREQUENCY(obj.CENTER_FREQUENCY >= freqStart & ...
%                 obj.CENTER_FREQUENCY <= freqEnd);
% 		end
% 	end
% 
% end
	
% function dataout = Octave(m, freStart, freEnd)
	% ==============================================================================
	% Provide list of octave frequencies and 1/3 octave frequencies.
	% ------------------------------------------------------------------------------
	% Input
	%	m		- m octave. 1 for octave and 3 for 1/3 octave
	%			- scalar
	%	freStart, freEnd
	%			- Start and end frequency in the range
	%			- scalar
	% ------------------------------------------------------------------------------
	% Output
	%	dataout.freCnt	- Center frequencies
	%					- Column
	%	dataout.freLow	- Low frequencies
	%					- Column
	%	dataout.freHi	- High frequencies
	%					- Column
	% ==============================================================================


	% if m == 1
		% freCnt = OCT_CNT;
    % elseif m == 3
		% freCnt = OCT_3_CNT;
	% end

	% The common joint factor due to m
	% factorM = sqrt(2 ^ (1 / m));

	% freHi = freCnt * factorM;
	% freLow = freCnt / factorM;

	% if nargin < 2
		% freStart = 25;
		% freEnd = 20000;
	% end

	% for i = 1 : length(freCnt)
		% if freStart <= freCnt(i)
			% break
		% end
	% end
	% for j = length(freCnt) : -1 : 1
		% if freEnd >= freCnt(j)
			% break
		% end
	% end
	% freCnt = freCnt(i : j);
	% freHi = freHi(i : j);
	% freLow = freLow(i : j);
	
	% dataout.freCnt = freCnt;
	% dataout.freHi = freHi;
	% dataout.freLow = freLow;
% end

