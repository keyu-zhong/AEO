clear
clc
close all;


N=120;  
Max_iter=200; 

all_cell_width = [231,308,385];
dim= 441;
Nt = 60;
rows = 21;
cols = 21;


for k=1:numel(all_cell_width)
    Cons = [];
    lb = 0.*ones(1,dim);
    ub = 1.*ones(1,dim);
    cell_width = all_cell_width(k);

    fobj = @(x) Fun(x, rows, cols, Nt, cell_width);

    Max_run = 1;
    X_init=zeros(N,dim);
    for kk=1:Max_run
        for m=1:N
           X_(m,:) = lb + rand(1,dim).*(ub - lb); 
        end
        X_T{kk}=X_;
    end
    for run = 1:Max_run
        X_init = X_T{run};
        %% AEO
        tic;
        [AEO_Best_Cost,AEO_Best_Power,AEO_Best_Efficiency,AEO_Pos,AEO_Cure] = AEO(X_init,N,Max_iter,lb,ub,dim,fobj,Cons,Nt);
        AEO_F1(k,run) = AEO_Best_Cost;
        AEO_F2(k,run) = AEO_Best_Power;
        AEO_F3(k,run) = AEO_Best_Efficiency;
        AEO_all_Pos{k,run} = AEO_Pos;
        AEO_all_cure{k,run} = AEO_Cure;
        AEO_time(k,run) = toc;
        
    end

end
save('./Results_3speeds_12_dirs/results');

