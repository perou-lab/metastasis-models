% `DWD_n5_Ensembl_log2_medcent_colstd.m`
% Use DWD to batch correct by dataset; adjusting all datasets to METABRIC
% reference

clear all

rng(123, 'twister')

addpath('BatchAdjust');

% Read in the training sets
trainp = readmatrix("METABRIC_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt");
trainp(:,1)=[];
trainn = readmatrix("G2034_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt");
trainn(:,1)=[];

% Generate direction vector between reference and second set
[dirvec, beta, dr] = DWD2XQ(trainp, trainn);

adjusted_trainp = trainp - dirvec * mean(trainp' * dirvec);
adjusted_trainn = trainn - dirvec * mean(trainn' * dirvec);

% Save direction vector
fileID = fopen("dirvec_METABRIC_G2034.txt","w");
fprintf(fileID, "%.12e\n", dirvec);
fclose(fileID);

% Read in new training set
trainn = readmatrix("G2603_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt");
trainn(:,1)=[];

% Generate direction vector between reference and second set
[dirvec, beta, dr] = DWD2XQ(trainp, trainn);

adjusted_trainp = trainp - dirvec * mean(trainp' * dirvec);
adjusted_trainn = trainn - dirvec * mean(trainn' * dirvec);

% Save direction vector
fileID = fopen("dirvec_METABRIC_G2603.txt","w");
fprintf(fileID, "%.12e\n", dirvec);
fclose(fileID); 

% Read in new training set
trainn = readmatrix("G12276_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt");
trainn(:,1)=[];

% Generate direction vector between reference and second set
[dirvec, beta, dr] = DWD2XQ(trainp, trainn);

adjusted_trainp = trainp - dirvec * mean(trainp' * dirvec);
adjusted_trainn = trainn - dirvec * mean(trainn' * dirvec);

% Save direction vector
fileID = fopen("dirvec_METABRIC_G12276.txt","w");
fprintf(fileID, "%.12e\n", dirvec);
fclose(fileID);

% Read in new training set
trainn = readmatrix("NKI295_Ensembl_clpsd_filt_imp_common_10814_log2_medcent_colsd.txt");
trainn(:,1)=[];

% Generate direction vector between reference and second set
[dirvec, beta, dr] = DWD2XQ(trainp, trainn);

adjusted_trainp = trainp - dirvec * mean(trainp' * dirvec);
adjusted_trainn = trainn - dirvec * mean(trainn' * dirvec);

% Save direction vector
fileID = fopen("dirvec_METABRIC_NKI295.txt","w");
fprintf(fileID, "%.12e\n", dirvec);
fclose(fileID);