clear
clc
close all;


N=120;  
Max_iter=200; 

cell_width = 77*3;
dim= 144;
Nt = 25;
rows = 12;
cols = 12;
%%%% 不同的单元格约束
Cons_L0 = [];  % L0
Cons_L1 = [121:144];  % L1
Cons_L2 = 61:84;   % L2
Cons_L3 = [11,12,23,24,35,36,47,48,59,60,71,72,83,84,95,96,107,108,119,120,131,132,143,144];   % L3
Cons_L4 = Cons_L3 - 5;   % L4
Cons_L5 = [41:44, 53:56, 65:68, 77:80, 89:92, 101:104];  % L5
Cons_L6 = [1,2,11,12,13,14,23,24,25,26,35,36,109,110,119,120,121,122,131,132,133,134,143,144];   % L6
Cons_L7 = 133:144;   % L7
Cons_L8 = 61:72;   % L8
Cons_L9 = [12,24,36,48,60,72,84,96,108,120,132,144];    % L9
Cons_L10 = Cons_L9 - 6;    % L10
Cons_L11 = [42,43,54,55,66,67,78,79,90,91,102,103]; % L11
Cons_L12 = [1,2,11,12,13,24,121,132,133,134,143,144];   % L12

for k=0:12
    
    lb = 0.*ones(1,dim);
    ub = 1.*ones(1,dim);
    eval(['Cons = Cons_L' num2str(k)]);
%     Cons = Cons_L7;
    ub(Cons) = 0;

    fobj = @(x) Fun(x, rows, cols, Nt, cell_width);

    X_init = zeros(N,dim);
    Max_run = 1;
     %% 初始化种群
     for kk=1:Max_run
        for m = 1:N
            X_(m,:) = lb + rand(1,dim).*(ub - lb);
        end
        X_T{kk} = X_;
     end
    for run = 1:Max_run
       X_init = X_T{run};
        %% AEO
        tic;
        [AEO_Best_Cost,AEO_Best_Power,AEO_Best_Efficiency,AEO_Pos,AEO_Cure] = AEO(X_init,N,Max_iter,lb,ub,dim,fobj,Cons,Nt);
        AEO_F1(k+1,run) = AEO_Best_Cost;
        AEO_F2(k+1,run) = AEO_Best_Power;
        AEO_F3(k+1,run) = AEO_Best_Efficiency;
        AEO_all_Pos{k+1,run} = AEO_Pos;
        AEO_all_cure{k+1,run} = AEO_Cure;
        AEO_time(k+1,run) = toc;
        
    end
    
end
save('./Results_1speeds_6_dirs/results');

