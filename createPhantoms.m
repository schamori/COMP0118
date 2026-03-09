function created_Image = createPhantoms(Model, TE, T2, S0, Nrow, Ncol, c_gt)
%%%% The function "Create_Phantoms" creates 256 * 256 phantom images without noise; Mean = 0, Sigma = 0
%%%% Inputs:-
%%%% Model:- "Exponential + Constant", "Exponential" or "Biexponential"
%%%% T2:- T2* value "One element"
%%%% Outputs:-
%%%% created_Image:- 3D image, the third dimension is for the intensity with Time to Echo "TE".
% ============================================================= %

% Pre-allocate the matrix for speed
created_Image = zeros(Nrow, Ncol, length(TE));

if strcmp(Model, 'exp+c')
    SI = S0 * exp(-TE / T2) + c_gt;                                         
    mkdir('exp_plus_c');
    Folder_Name = 'exp_plus_c\';
elseif  strcmp(Model, 'exp')
    SI = S0 * exp(-TE / T2);                 
    mkdir('exp');
    Folder_Name = 'exp\';
elseif strcmp(Model, 'exp+exp') % Fixed syntax error here
    SI = 0.9 * S0 * exp(-TE / T2) + 0.1 * S0 * exp(-TE / 200);    
    mkdir('exp_exp');
    Folder_Name = 'exp_exp\';
else
    error('Unknown model type'); % Added error handling
end

for i = 1:length(TE)
    created_Image(:,:,i) = SI(i) .* ones(Nrow, Ncol);
    
    % Example of how to use the Folder_Name to save data:
    % save([Folder_Name, 'image_', num2str(i), '.mat'], 'created_Image');
end
end