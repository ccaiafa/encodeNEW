function [] = split_shells(dwiFile,bvecsFile,bvalsFile,b0_max,bvals_rounded,shell)
% normalizes the bvals and splits the bvecs

[path,~,~] = fileparts(fullpath(dwiFile));

% Parameters used for normalization
bvals = dlmread(bvalsFile);
bvecs = dlmread(bvecsFile);
dwi = niftiRead(dwiFile);

% Round the numbers to the closest thousand 
[bvals_unique, ~, bvals_uindex] = unique(bvals);
bvals_unique(bvals_unique <= b0_max) = 0;
bvals_unique = round(bvals_unique./bvals_rounded)*bvals_rounded;
bvals_round = bvals_unique( bvals_uindex );

index = (bvals_round == shell);
index0 = (bvals_round == 0);
all_index = or(index, index0);
assertEqual(sum(all_index), sum(index0) + sum(index));

%write files
dlmwrite(strcat(path,'/dwi_',num2str(shell),'.bvals'), bvals_round(all_index), 'delimiter',' ');
dlmwrite(strcat(path,'/dwi_',num2str(shell),'.bvecs'), bvecs(:,all_index), 'delimiter', ' ');
dwi_oneshell = dwi;
dwi_oneshell.fname = strcat(path,'/dwi_',num2str(shell),'.nii.gz');
dwi_oneshell.data = dwi.data(:,:,:,all_index);
dwi_oneshell.dim(4) = size(dwi_oneshell.data,4);
niftiWrite(dwi_oneshell);