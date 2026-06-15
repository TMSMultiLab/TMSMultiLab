%% LOAD A MATRIX OF COLORCODES TO USE IN MATLAB____________________________
% Copyright (c) 2014, Diana Di Leonardo All rights reserved.
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESSINTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

colorblind =...
[0.0000    0.0000    1.0000;
 1.0000    0.0000    0.0000;
 1.0000    1.0000    0.0000;
 0.6602    0.6602    0.6602;
 0.0000    0.0000    0.0000;
 1.0000    0.6445    0.0000;
 1.0000    0.0000    1.0000;
 0.0000    0.5000    0.5000;
 0.0000    0.0000    0.5430;
 0.0000    0.3906    0.0000;
 0.0000    1.0000    1.0000;
 0.5977    0.1953    0.7969];
colornames = {'blue';'red';'yellow';'darkgray';'black';'orange';'magenta';'teal';'darkblue';'darkgreen';'cyan';'darkorchid'};

colorblind7 =... % for 7-levels, approximately from reds to blues
[0.75,0.20,0.33;
 1.00,0.50,0.25;
 1.00,0.75,0.17;
 0.50,0.50,0.50;
 0.50,0.50,0.83;
 0.00,0.75,0.67;
 0.00,0.50,1.00;];
% red			rgb(191, 51, 84)	#bf3354
% dark orange		rgb(255, 128, 51)	#ff8033
% light orange		rgb(255, 191, 43)	#ffbf2b
% mid grey		rgb(128, 128, 128)	#808080
% mauve			rgb(128, 128, 212)	#8080d4 
% teal			rgb(0, 191, 171)	#00bfab
% blue			rgb(0, 128, 255)	#0080ff