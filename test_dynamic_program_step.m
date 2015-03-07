%% This function generates a system with different terminal functions and constraints but
% with same size. The constraint are preconditioned accodingly.
clear all;
close all;
clc;
Nm=5; % Number of masses
T_sampling=0.5;
ops_masses=struct('Ts',T_sampling,'xmin', ...
    -4*ones(2*Nm,1), 'xmax', 4*ones(2*Nm,1), 'umin', -2*ones(Nm-1,1),'umax',...
    2*ones(Nm-1,1),'b', 0.1*ones(Nm+1,1));
ops_system.nx=2*Nm;
ops_system.nu=Nm-1;
ops_system.sys_uncert=1;
ops_system.ops_masses=ops_masses;
%sys_no_precond=system_masses(Nm,ops_masses);
predict_horz=8;%prediction horizon
Test_points=10;
x_rand=4*rand(ops_system.nx,Test_points)-2;
result.u=zeros(ops_system.nu,Test_points);
time_gpad=cell(Test_points,1);
time_gurobi=cell(Test_points,1);
U_max=zeros(2,Test_points);
U_min=zeros(2,Test_points);
dual_gap=zeros(2,Test_points);
%% Generation of tree
scenario_size=[5 3 2];
%ops_system.Np=various_predict_horz(no_of_pred);
ops.N=predict_horz; %step 2: argmin of the lagrangian using dynamic programming
ops.brch_ftr=ones(ops.N,1);
ops.brch_ftr(1:size(scenario_size,2))=scenario_size(1:size(scenario_size,2));
Ns=prod(ops.brch_ftr);
ops.nx=ops_system.nx;
ops.prob=cell(ops.N,1);
for i=1:ops.N;
    if(i<=size(scenario_size,2))
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
ops_GPAD.x0=x_rand(:,1);
ops_GPAD.primal_inf=1e-3;
ops_GPAD.dual_gap=10e-3;
ops_GPAD.alpha=1/calculate_Lipschitz(sys,V,Tree);

Ns=length(Tree.leaves);% total scenarios in the tree
Nd=length(Tree.stage);%toal nodes in the tree

for kk=1:Test_points
    
    % Initalizing the dual varibables
    Y.y=0.1*rand(Nd-Ns,size(sys.F{1},1));
    for i=1:Ns
        Y.yt{i,:}=0.1*rand(1,size(sys.Ft{i,1},1));
    end
    
    [Z,Q]=GPAD_dynamic_multiple(sys,Ptree,Tree,Y,ops_GPAD.x0);
    
    yalmip_dp=yalmip_dynamic_step(sys,V,Tree,Y);%yalmip variable
    Z_yalmip=yalmip_dp{ops_GPAD.x0};
    
    prm_cst=0;
    for i=1:Tree.ancestor(Tree.leaves(end))
        prm_cst=prm_cst+Tree.prob(i)*(Z.X(:,i)'*V.Q*Z.X(:,i)+Z.U(:,i)'*V.R*Z.U(:,i))...
            +Y.y(i,:)*(sys.F{i}*Z.X(:,i)+sys.G{i}*Z.U(:,i)-sys.g{i});
    end
    
    for i=1:length(Tree.leaves)
        prm_cst=prm_cst+Tree.prob(Tree.leaves(i))*(Z.X(:,Tree.leaves(i))'*V.Vf{i}*...
            Z.X(:,Tree.leaves(i)))+Y.yt{i}*(sys.Ft{i}*Z.X(:,Tree.leaves(i))-sys.gt{i});
    end
    [Z_yalmip{1,3} prm_cst]
    Z_yalmip{1,3}-prm_cst
    %[max(max(Z_yalmip{1,2}(:,:)-Z.U(:,:))) min(min(Z_yalmip{1,2}(:,:)-Z.U(:,:)))]
    %[max(max(Z_yalmip{1,1}(:,:)-Z.X(:,:))) min(min(Z_yalmip{1,1}(:,:)-Z.X(:,:)))]
    
    prm_cst1=0;
    Z1.X=Z_yalmip{1,1}(:,:);
    Z1.U=Z_yalmip{1,2}(:,:);
    for i=1:Tree.ancestor(Tree.leaves(end))
        prm_cst1=prm_cst1+Tree.prob(i)*(Z1.X(:,i)'*V.Q*Z1.X(:,i)+Z1.U(:,i)'*V.R*Z1.U(:,i))...
            +Y.y(i,:)*(sys.F{i}*Z1.X(:,i)+sys.G{i}*Z1.U(:,i)-sys.g{i});
    end
    
    for i=1:length(Tree.leaves)
        prm_cst1=prm_cst1+Tree.prob(Tree.leaves(i))*(Z1.X(:,Tree.leaves(i))'*V.Vf{i}*...
            Z1.X(:,Tree.leaves(i)))+Y.yt{i}*(sys.Ft{i}*Z1.X(:,Tree.leaves(i))-sys.gt{i});
    end
    prm_cst-prm_cst1
   % x_rand=4*rand(ops_system.nx,1)-2;
end