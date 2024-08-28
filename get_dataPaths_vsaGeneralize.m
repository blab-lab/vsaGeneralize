function [dataPaths] = get_dataPaths_vsaGeneralize(bPilot)
if nargin < 1
    bPilot = 0;
end

if bPilot
    svec = {'test_1011_JE',...
            'test_1016_SS',...
            'test_1030_SA',...
            'test_1101_HS',...
            'test_1104_MK',...
            };
else
    %don't include 252 because of crash
    svec = [146 189 235 253 255 257 258 266 268 269 277 280 285 289 291 292 294 295 296 298 299];
end

dataPaths = get_acoustLoadPaths('vsaGeneralize',svec);
