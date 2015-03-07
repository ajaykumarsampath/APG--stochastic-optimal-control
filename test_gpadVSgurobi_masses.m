%%
% This function generates a system with different terminal functions and constraints but
% with same size. The constraint are preconditioned accodingly.
% We solve the method using different methods. First formulated using
% 1) Gurobi-IP 2) Gurobi-AS directly

clear all;
close all;
clear model;
clc;
Nm=5; % Number of masses
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
Test_points=1;
x_rand=4*rand(ops_system.nx,Test_points)-2;
time_gpad=cell(Test_points,1);
U_max=zeros(2,Test_points);
U_min=zeros(2,Test_points);
%dual_gap=zeros(2,Test_points);
test_cuda=0;
%% Yalmip details

details_solvers.IP=1;
details_solvers.method=details_solvers.IP*1;
time_solver=zeros(Test_points,details_solvers.method+2);
result.u=cell(details_solvers.method+1,1);

params.outputflag=0;
%% Generation of tree
scenario_size=[2 1 1];
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
        %Cost functioin
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
            %consitraint in the horizon
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
        
        [sys,Hessian_iapp]=calculate_diffnt_precondition_matrix(sys_no_precond,V,Tree...
            ,struct('use_cell',1,'use_hessian',0));
        tic;
        Ptree=GPAD_dynamic_formul_precond_multip(sys,V,Tree);
        toc
        
        ops_GPAD.steps=200;
        ops_GPAD.primal_inf=1e-3;
        ops_GPAD.dual_gap=10e-3;
        ops_GPAD.alpha=1/calculate_Lipschitz(sys,V,Tree);
        %}
        if(details_solvers.method)
            
            %GUROBI-IP
            [yalmip_tree{1}, yalmip_const_time{1}]=yalmip_primal_multiple(sys_no_precond,V,Tree);%yalmip variable
            
        end
        
        model=gurobi_solve(sys_no_precond,V,Tree);
        %}
        max_size=zeros(Test_points,length(Tree.stage));
        
        for kk=1:Test_points
            
            ops_GPAD.x0=x_rand(:,kk);
            if(details_solvers.method)
                tic;
                Z_yalmip=yalmip_tree{1}{ops_GPAD.x0};
                time_solver(3,kk)=toc;
                result.u{1}(:,kk)=Z_yalmip{1,2}(:,1);
            end
            %}
            model.rhs(end-sys_no_precond.nx+1:end)=ops_GPAD.x0;
            %gurobi_write(model, 'sc_qp.lp');
            tic
            results=gurobi(model,params);
            time_solver(1,kk)=toc;
            time_solver(2,kk)=results.runtime;
            nvar=sys_no_precond.nx+sys_no_precond.nu;
            non_leave_nodes=size(Tree.children,1);
            for jj=1:length(Tree.prob)
                if(jj<=non_leave_nodes)
                    Zgurobi{1,1}(:,jj)=results.x((jj-1)*nvar+1:(jj-1)*nvar+sys_no_precond.nx,1);
                    Zgurobi{1,2}(:,jj)=results.x((jj-1)*nvar+sys_no_precond.nx+1:jj*nvar,1);
                else
                    Zgurobi{1,1}(:,jj)=results.x(non_leave_nodes*nvar+(jj-non_leave_nodes-1)...
                        *sys_no_precond.nx+1:non_leave_nodes*nvar+(jj-non_leave_nodes)*sys_no_precond.nx);
                end
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
%transfer_multiple_data


figure(1)
plot(time_solver(1,:));
hold all;
plot(time_solver(2,:));
plot(time_solver(3,:));

mean(time_solver(1,:))
mean(time_solver(2,:))

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