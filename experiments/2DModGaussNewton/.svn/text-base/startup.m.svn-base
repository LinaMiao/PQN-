% read SLIM_ROOT environment
slimroot = getenv('SLIM_ROOT');
if length(slimroot)<1
	fprintf('FATAL ERROR: SLIM_ROOT environment does not exist\n');
	fprintf('\t source appropriate environment.(sh/csh) in installation directory\n');
	break;
end

% start
fprintf('Loading SLIM Toolboxes from\n\t%s\n',slimroot);

% SPOT (slim updates)
addpath([slimroot '/tools/utilities/SPOT-SLIM']);
addpath([slimroot '/tools/utilities/SPOT-SLIM/tests']);
try
   addpath(fullfile(spot.path,'tests','xunit'));
catch ME
   error('Can''t find xunit toolbox.');
end
% pSPOT
addpath([slimroot '/tools/utilities/pSPOT']);
% SPGL1-SLIM (parallel)
addpath([slimroot '/tools/solvers/SPGL1-SLIM']);
% CurveLab
addpath([slimroot '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_usfft_matlab']);
addpath([slimroot '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_matlab']);
addpath([slimroot '/tools/transforms/CurveLab-2.1.2-SLIM/fdct_wrapping_cpp/mex']);
addpath([slimroot '/tools/transforms/CurveLab-2.1.2-SLIM/fdct3d/mex']);
addpath([slimroot '/tools/transforms/CurveLab-2.1.2-SLIM/mecv']);
% REPSI
addpath([slimroot '/tools/algorithms/REPSI']);
% 2DFreqModeling
addpath([slimroot '/tools/algorithms/2DFreqModeling']);
% Miscellaneous
addpath([slimroot '/tools/operators/misc']);
addpath([slimroot '/tools/functions/misc']);

% MADAGASCAR
addpath([slimroot '/external/lib']);

% done
fprintf('Done loading SLIM Toolboxes\n');
%path
addpath(genpath(pwd));