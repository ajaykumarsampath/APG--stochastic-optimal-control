%%
% This function generates a system with different terminal functions and constraints but
% with same size. The constraint are preconditioned accodingly.
% We solve the method using different methods. First formulated using
% 1) Gurobi-IP 2) Gurobi-AS 3) qpOASES 4)QPC-IP 5)QPC-AS

clear all;
close all;
clc;
Nm=2; % Number of masses
T_sampling=0.5;
ops_masses=struct('Ts',T_sampling,'xmin', ...
    -4*ones(2*Nm,1), 'xmax', 4*ones(2*Nm,1), 'umin', -2*ones(Nm-1,1),'umax',...
    2*ones(Nm-1,1),'b', 0.1*ones(Nm+1,1));
ops_system.nx=2*Nm;
ops_system.nu=Nm-1;
ops_system.sys_uncert=0;
ops_system.ops_masses=ops_masses;
%sys_no_precond=system_masses(Nm,ops_masses);
various_predict_horz=10;%prediction horizon
Test_points=100;
x_rand=4*rand(ops_system.nx,Test_points)-2;
time_gpad=cell(Test_points,1);
U_max=zeros(2,Test_points);
U_min=zeros(2,Test_points);
%dual_gap=zeros(2,Test_points);
test_cuda=0;
%% Yalmip details

details_solvers.IP=1;
details_solvers.AS=0; % Active-Set details
details_solvers.QPC=0; % QPC solver use
details_solvers.qpoases=0;% qpoases solver
details_solvers.method=details_solvers.IP*1+(2+details_solvers.QPC)*details_solvers.AS+...
    details_solvers.QPC;
time_solver=cell(Test_points,details_solvers.method+1);
result.u=cell(details_solvers.method+1,1);
%% Generation of tree
scenario_size=[2 2 1];
for N_prb_steps=3:length(scenario_size)
    for no_of_pred=1:length(various_predict_horz)
        %ops_system.Np=various_predict_horz(no_of_pred);
        ops.N=various_predict_horz(no_of_pred); %step 2: argmin of the lagrangian using dynamic programming
        ops.brch_ftr=ones(ops.N,1);
        ops.brch_ftr(1:N_prb_steps)=scenario_size(1:N_prb_steps);
        Ns=prod(ops.brch_ftr);
        ops.nx=ops_system.nx;
        ops.prob=cell(ops.N,1);
        for i=1:ops.N;
            if(i<=N_prb_steps)
                pd=rand(1,ops.brch_ftr(i));
                if(i==1)
                    ops.prob{i,1}=pd/sum(pd);
                    pm=1;
                else
                    pm=pm*scenario_size(i-1);
                    ops.prob{i,1}=kron(ones(pm,1),pd/sum(pd));
                end
            else
                ops.prob{i,1}=ones(1,Ns);
            end
        end
        tic
        [sys_no_precond,Tree]=tree_generation_multiple(ops_system,ops);
        time.tree_formation=toc;
        sys_no_precond.nx=ops_system.nx;
        sys_no_precond.nu=ops_system.nu;
        sys_no_precond.Np=ops.N;
        SI=scenario_index(Tree);%calculation of the scenario index.
        %%
        %Cost function
        V.Q=eye(sys_no_precond.nx);
        V.R=eye(sys_no_precond.nu);
        %%terminal constraints
        sys_no_precond.Ft=cell(Ns,1);
        sys_no_precond.gt=cell(Ns,1);
        V.Vf=cell(Ns,1);
        sys_no_precond.trm_size=(2*sys_no_precond.nx)*ones(Ns,1);
        %r=rand(Ns,1);
        r=ones(Ns,1);
        for i=1:Ns
            %constraint in the horizon
            sys_no_precond.Ft{i}=[eye(sys_no_precond.nx);-eye(sys_no_precond.nx)];
            sys_no_precond.gt{i}=(3+0.1*rand(1))*ones(2*sys_no_precond.nx,1);
            nt=size(sys_no_precond.Ft{i},1);
            P=Polyhedron('A',sys_no_precond.Ft{i},'b',sys_no_precond.gt{i});
            if(isempty(P))
                error('Polyhedron is empty');
            end
            V.Vf{i}=dare(sys_no_precond.A{1},sys_no_precond.B{1},r(i)*V.Q,r(i)*V.R);
            %V.Vf{i}=dare(sys_no_precond.A,sys_no_precond.B,r(i)*V.Q,r(i)*V.R);
        end
        %% preconditioning the system and solve the system using dgpad.
        [sys,Hessian_app]=calculate_diffnt_precondition_matrix(sys_no_precond,V,Tree...
            ,struct('use_cell',1,'use_hessian',0));
        tic;
        Ptree=GPAD_dynamic_formul_precond_multip(sys,V,Tree);
        toc
        ops_GPAD.steps=200;
        ops_GPAD.primal_inf=1e-3;
        ops_GPAD.dual_gap=10e-3;
        ops_GPAD.alpha=1/calculate_Lipschitz(sys,V,Tree);
        
        if(details_solvers.method)
            
            %GUROBI-IP
            [yalmip_tree{1}, yalmip_const_time{1}]=yalmip_primal_multiple(sys_no_precond,V,Tree);%yalmip variable
            
            if(details_solvers.AS)
                %GUROBI-AS
                tree_solver_settings= sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
                tree_solver_settings.gurobi.method=0;
                [yalmip_tree{details_solvers.AS*2}, yalmip_const_time{2}]=yalmip_primal_multiple(sys_no_precond,V,Tree...
                    ,tree_solver_settings);%yalmip variable
                nvar*child+K*sys.nx
                if(details_solvers.qpoases)
                    %qpOASES
                    tree_solver_settings= sdpsettings('solver','qpoases','verbose',0,'cachesolvers',1);
                    tree_solver_settings.qpoases=qpOASES_options('mpc');
                    tree_solver_settings.qpoases.printLevel=1;
                    p=details_solvers.IP+details_solvers.AS+details_solvers.qpoases;
                    [yalmip_tree{p}, yalmip_const_time{3}]=yalmip_primal_multiple(sys_no_precond,V,Tree,...
                        tree_solver_settings);%yalmip variable
                end
            end
            %QCP
            if(details_solvers.QPC)
                
                %IP solver
                tree_solver_settings= sdpsettings('solver','qpip','verbose',0,'cachesolvers',1);
                p=details_solvers.AS*2+details_solvers.IP+details_solvers.QPC;
                [yalmip_tree{p}, yalmip_const_time{2}]=yalmip_primal_multiple(sys_no_precond,V,Tree...
                    ,tree_solver_settings);%yalmip variable
                
                if(details_solvers.AS)
                    %AS solver
                    p=details_solvers.AS*3+2;
                    tree_solver_settings= sdpsettings('solver','qpas','verbose',0,'cachesolvers',1);
                    [yalmip_tree{p}, yalmip_const_time{2}]=yalmip_primal_multiple(sys_no_precond,V,Tree...
                        ,tree_solver_settings);%yalmip variable
                end
            end
        end
        %}
        max_size=zeros(Test_points,length(Tree.stage));
        
        for kk=1:Test_points
            
            ops_GPAD.x0=x_rand(:,kk);
            for jj=1:details_solvers.method
                tic;
                Z_yalmip=yalmip_tree{jj}{ops_GPAD.x0};
                time_solver{kk,jj}.yalmip_time=toc;
                result.u{jj}(:,kk)=Z_yalmip{1,2}(:,1);
                time_solver{kk,jj}.yalmip_cost=Z_yalmip{1,3};
            end
            
            %{
            %GPAD
            [Z_gpad_pre,Y_gpad_pre,time_gpad{kk}]=GPAD_multiple(sys,Ptree,Tree,V,ops_GPAD);
            if(~isfield(time_gpad{kk},'iterate'))
                time_gpad{kk}.iterate=ops_GPAD.steps;
            end
            %max(max(Z_yalmip{1,2}(:,:)-Z_gpad_pre.U(:,:)))
            U_max(2,kk)=max(max(Z_gpad_pre.U-Z_yalmip{1,2}));
            U_min(2,kk)=min(min(Z_gpad_pre.U-Z_yalmip{1,2}));
            %dual_gap(2,kk)=time_gpad{kk}.dual_gap;
            %}
        end
    end
end

%% testing the dual gradient calculation
if(test_cuda)
    x_dgpad=rand(sys.nx,1);
    Ns=length(Tree.leaves);% total scenarios in the tree
    Nd=length(Tree.stage);
    W.y=rand(Nd-Ns,size(sys.F{1},1));
    for i=1:Ns
        W.yt{i,:}=rand(1,size(sys.Ft{i,1},1));
    end
    Zdgpad=GPAD_dynamic_multiple(sys,Ptree,Tree,W,x_dgpad);
    
    yalmip_dp=yalmip_dynamic_step(sys,V,Tree,W);
    Zyalmip_dp=yalmip_dp{x_dgpad};
    assert(norm(Zyalmip_dp{1,1}-Zdgpad.X,inf)<1e-4);
    assert(norm(Zyalmip_dp{1,2}-Zdgpad.U,inf)<1e-4);
end
%%
transfer_multiple_data


iterates=zeros(1,Test_points);
time_operation=zeros(6,Test_points);
dual_gap=zeros(1,Test_points);
prim_cost=zeros(6,Test_points);
cst_cmp=zeros(1,Test_points);
for i=1:Test_points
    for k=1:details_solvers.method
        time_operation(k+1,i)=time_solver{i,k}.yalmip_time;
        prim_cost(k+1,i)=time_solver{i,k}.yalmip_cost;
    end
end

if(details_solvers.method)
    plot(1);
    for k=1:details_solvers.method
        plot(time_operation(k+1,2:end));
        hold all;
    end
    legend('Gurobi-IP','Gurobi(PAC)','qpOASES');
end
mean(time_operation(2,2:end))
max(time_operation(2,2:end))
%{
for i=1:Test_points
    iterates(1,i)=time_gpad{i,1}.iterate;
    time_operation(1,i)=time_gpad{i,1}.gpad_solve;
    time_operation(2,i)=time_gurobi{i}.yalmip_time;
    dual_gap(1,i)=time_gpad{i,1}.dual_gap(end);
    cst_cmp(1,i)=time_gpad{i}.prm_cst(end)-time_gurobi{i}.yalmip_cost;
    prim_cost(1,i)=time_gpad{i}.prm_cst(end);
    prim_cost(2,i)=time_gurobi{i}.yalmip_cost;
end
mean(time_operation')
figure(1)
subplot(211)
plot(iterates)
title('Iterates')
subplot(212)
plot(time_operation(1,:))
title('Gpad time')
figure(2)
subplot(211)
plot(cst_cmp)
title('difference of the cost');
subplot(212)
plot(prim_cost(1,:))
hold all;
plot(prim_cost(2,:))
title('Actual cost')

%}
%mean(time_operation')
%max(time_operation')
%mean(iterates')
%max(iterates')