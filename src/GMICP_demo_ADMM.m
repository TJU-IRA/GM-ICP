%% 基于ADMM方法的GM-ICP demo  by yugeng.huang
clc;clear;close all;

    source_file = 'Bunny\sourcefile.ply';
    target_file = 'Bunny\targetfile.ply';
    gt = load('Bunny\gt.txt');

Ansfile =  'p0.3.txt';

alpha1 = 2.0; %% 权重参数 默认2
PP_norm = 1.0; %% 阶次 默认1.0 性能最优

writeAsns = false;
Tini_gt = gt;

source = pcread(source_file);
target = pcread(target_file);

SP = double(source.Location');    %source points
SN = double(source.Normal');      %source normals
TP = double(target.Location');    %target points
TN = double(target.Normal');      %target normals
figure
pcshow(SP', 'r'), 
hold on
pcshow(TP','g')
%% 点云归一化
scaleS = norm(max(SP,[],2)-min(SP,[],2));
scaleT = norm(max(TP,[],2)-min(TP,[],2));
scale = max(scaleS,scaleT);

SP = SP/scale;
TP = TP/scale;
    %%
%%
t1 = clock;
[T0, count]= GM_ICP_ADMM(SP, TP, SN, TN, alpha1,PP_norm);
t2 = clock;
time = etime(t2,t1);
trans = T0(1:3,4);
trans = trans*scale;
T0(1:3,4)=trans;
SP = double(source.Location');
P1 = T0(1:3,1:3)*SP+repmat(T0(1:3,4),1,size(SP,2));
P2 = Tini_gt(1:3,1:3)*SP+repmat(Tini_gt(1:3,4),1,size(SP,2));
rmse = sqrt(sum(sum((P1-P2).^2))/size(SP,2))

TP = double(target.Location'); 
figure
pcshow(P1', 'r'), hold on
pcshow(TP','g')
GlobalCompare2(T0, Tini_gt);

if writeAsns
    fid = fopen(Ansfile,'w');
    fprintf(fid,[repmat('%5.6f\t', 1, size(T0,2)), '\n'], T0');
    fclose(fid);
end


