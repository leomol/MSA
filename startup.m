function startup()
    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(fullfile(root, 'common'));
end