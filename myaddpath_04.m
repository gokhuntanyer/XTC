%% test_addpath_MyFcns
% add the directory as given below, so you can 
% call all the functions in the directory MyFcns
% S. G. Tanyer
% 180211
%
% Include the following lines at the beginning of your code. 
% directories;
% MATLAB-MyFcns
% DATA-180307-TXdata
%
% mydir1 = cd;
% cd ../
% cd ../
% mydir2 = cd;
% cd MATLAB-MyFcns
% mydir3 = cd;
% addpath (mydir1, mydir2, mydir3);
% fprintf('____________________________________\n');
% fprintf('Directory included in the path: \n');
% fprintf('>%s\n' , mydir3);
% fprintf('Current Directory is: \n');
% fprintf('>%s\n' , mydir1);
clc;
mydir1 = cd;
%
cd ../
cd  DATA-180331-HRTFs
mydir2 = cd;
%
cd ../
cd  DATA-180326-TX-Recdata
%cd DATA-180307-TXdata
mydir3 = cd;
%
cd ../
cd ../
cd MATLAB-MyFcns
%
mydir4 = cd;
addpath (mydir1, mydir2, mydir3, mydir4);
cd(mydir1);
%
fprintf('____________________________________\n');
fprintf('Directory included in the path: \n');
fprintf('>%s\n' , mydir2);
fprintf('>%s\n\n' , mydir3);
fprintf('Current Directory is: \n');
fprintf('>%s' , mydir1);
%
mytic;
%% THE END  THE END  THE END  THE END  THE END  THE END  THE END 