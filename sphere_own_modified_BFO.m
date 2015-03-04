% Rosenbrock minimization using BFO %
% by Kalyan Sourav Dash  %

clear all;
close all;
clc;

% Bacteria Foraging Optimization %

N=2;
n=1;
alpha=0.5;
beta=0.08;

% ------- initialisation ----------%
Ne=20;
Nr=20;
Nc=20;
Np=20;
Ns=10;
D=5;
C=0.01;
C_array=zeros(1,N);
C_array(1)=C;
Ped=0.9; % elimination dispersion probability

x=(rand(Np,D)-0.5)*60; % x lies in [-30 30]

J=zeros(Np,1);

for k=1:Np
    for i=1:D
        
    J(k)=sum(x(k,i)^2); % initial fitness calculation
    
    end
end
Jlast=J;

for w=1:N
    for l=1:Ne
        for k=1:Nr
            Jchem=J;
            for j=1:Nc
                % Chemotaxis Loop %

                for i=1:Np
                    del=(rand(1,D)-0.5)*2;
                    x(i,:)=x(i,:)+(C/sqrt(del*del'))*del;
                    for d=1:D
                        J(i)=sum(x(i,d)^2);
                    end

                    for m=1:Ns
                        if J(i)<Jlast(i)
                            Jlast(i)=J(i);
                            x(i,:)=x(i,:)+C*(del/sqrt(del*del'));
                            for d=1:D
                                J(i)=sum(x(i,d)^2);    
                            end
                        else
                            del=(rand(1,D)-0.5)*2;
                            x(i,:)=x(i,:)+C*(del/sqrt(del*del'));
                            for d=1:D
                                J(i)=sum(x(i,d)^2);    
                            end
                        end   
                    end

                end

                Jchem=[Jchem J];
            end  % End of Chemotaxis %


                for i=1:Np
                    Jhealth(i)=sum(Jchem(i,:)); % sum of cost function of all chemotactic loops for a given k & l
                end
                [Jhealth1,I]=sort(Jhealth,'ascend');
                x=[x(I(1:Np/2),:);x(I(1:Np/2),:)];
                J=[J(I(1:Np/2),:);J(I(1:Np/2),:)];
                xmin=x(I(1),:);
        end
        Jmin(l)=min(J);
        % random elimination dispersion

      for i=1:Np
            r=rand;
            if r>=Ped
                x(i,:)=(rand(1,D)-0.5);
                for d=1:D
                        J(i)=sum(x(i,d)^2);
                end
            end
        end

    end
    if w<=n
        C=0.01;
        C=C_array(w);
    else
        C_array(w)=C_array(w-1)*alpha;
        C=C_array(w);
    end
end
% plot(Jmin);
J_value_array=sort(Jmin,'ascend');
J_value_array(1)
                    