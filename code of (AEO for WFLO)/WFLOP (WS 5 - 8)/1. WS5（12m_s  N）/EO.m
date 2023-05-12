%_________________________________________________________________________________
%  Equilibrium Optimizer source code (Developed in MATLAB R2015a)
%
%  programming: Afshin Faramarzi & Seyedali Mirjalili
%
%  e-Mail: afaramar@hawk.iit.edu, afshin.faramarzi@gmail.com
%
%  paper:
%  A. Faramarzi, M. Heidarinejad, B. Stephens, S. Mirjalili, 
%  Equilibrium optimizer: A novel optimization algorithm
%  Knowledge-Based Systems
%  DOI: https://doi.org/10.1016/j.knosys.2019.105190
%____________________________________________________________________________________

function [Best_F1,Best_F2,Ceq1_fit,Ceq1,Cure] = EO(C,Particles_no,Max_iter,lb,ub,dim,fobj,Cons,Nt)

Ceq1=zeros(1,dim);   Ceq1_fit=inf; 
Ceq2=zeros(1,dim);   Ceq2_fit=inf; 
Ceq3=zeros(1,dim);   Ceq3_fit=inf; 
Ceq4=zeros(1,dim);   Ceq4_fit=inf;

% C=initialization_EO(Particles_no,dim,ub,lb);

Iter=0; V=1;

a1=2;
a2=1;
GP=0.5;

count = 1;
objective = inf.*ones(1,3);

while Iter<Max_iter
   
      for i=1:size(C,1)  
        
        Flag4ub=C(i,:)>ub;
        Flag4lb=C(i,:)<lb;
        C(i,:)=(C(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %%% ��ɢ��0��1���޸�����
        C(i,:) = round(C(i,:));
        C(i,:) = Repair(C(i,:), Nt, Cons);
        
        [f1,f2,fitness(i),power_order(i,:)] = fobj(C(i,:));

      
        if fitness(i)<Ceq1_fit 
              Ceq1_fit=fitness(i);  Ceq1=C(i,:);Best_F1 = f1;Best_F2 = f2;
        elseif fitness(i)>Ceq1_fit && fitness(i)<Ceq2_fit  
              Ceq2_fit=fitness(i);  Ceq2=C(i,:);              
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)<Ceq3_fit
              Ceq3_fit=fitness(i);  Ceq3=C(i,:);
        elseif fitness(i)>Ceq1_fit && fitness(i)>Ceq2_fit && fitness(i)>Ceq3_fit && fitness(i)<Ceq4_fit
              Ceq4_fit=fitness(i);  Ceq4=C(i,:);
                         
        end
      end
      
%---------------- Memory saving-------------------   
      if Iter==0
        fit_old=fitness;  C_old=C;
      end
    
     for i=1:Particles_no
         if fit_old(i)<fitness(i)
             fitness(i)=fit_old(i); C(i,:)=C_old(i,:);
         end
     end

    C_old=C;  fit_old=fitness;
%-------------------------------------------------
       
Ceq_ave=(Ceq1+Ceq2+Ceq3+Ceq4)/4;                              % averaged candidate 
C_pool=[Ceq1; Ceq2; Ceq3; Ceq4; Ceq_ave];                     % Equilibrium pool

 
 t=(1-Iter/Max_iter)^(a2*Iter/Max_iter);                      % Eq (9)

 
    for i=1:Particles_no
           lambda=rand(1,dim);                                % lambda in Eq(11)
           r=rand(1,dim);                                     % r in Eq(11)  
           Ceq=C_pool(randi(size(C_pool,1)),:);               % random selection of one candidate from the pool
           F=a1*sign(r-0.5).*(exp(-lambda.*t)-1);             % Eq(11)
           r1=rand(); r2=rand();                              % r1 and r2 in Eq(15)
           GCP=0.5*r1*ones(1,dim)*(r2>=GP);                   % Eq(15)
           G0=GCP.*(Ceq-lambda.*C(i,:));                      % Eq(14)
           G=G0.*F;                                           % Eq(13)
           C(i,:)=Ceq+(C(i,:)-Ceq).*F+(G./lambda*V).*(1-F);   % Eq(16)                                                             
    end
 
       Iter=Iter+1;  
       disp(['Iteration ' num2str(Iter) ': Best Cost = ' num2str(Ceq1_fit)]);
       Convergence_curve(1,Iter)=Ceq1_fit; 
       cure_cost_per(1,Iter) = Best_F1;
       cure_Farmpower(1,Iter) = Best_F2;
end
Cure = [cure_cost_per;cure_Farmpower;Convergence_curve];

display(['The best Cost_per found by EO is : ', num2str(Best_F1,10)]);
display(['The best Farmpower found by EO is : ', num2str(Best_F2,10)]);
display(['The best Efficiency found by EO is : ', num2str(Ceq1_fit,10)]);
fprintf('--------------------------------------\n');
end



function [Cin,domain]=initialization_EO(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Cin=rand(SearchAgents_no,dim).*(ub-lb)+lb;
    domain=ones(1,dim)*(ub-lb);
end


% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Cin(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
    domain=ones(1,dim).*(ub-lb);
end
end
