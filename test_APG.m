%%
% This function generates a system with different terminal functions and constraints but
% with same size. The constraint are preconditioned accodingly.

clear all;
close all;
clc;
Nm=2; % Number of masses
T_sampling=0.5; % options of the system.
ops_masses=struct('Ts',T_sampling,'xmin', ...
    -4*ones(2*Nm,1), 'xmax', 4*ones(2*Nm,1), 'umin', -2*ones(Nm-1,1),'umax',...
    2*ones(Nm-1,1),'b', 0.1*ones(Nm+1,1));

ops_system.nx=2*Nm;
ops_system.nu=Nm-1;
ops_system.sys_uncert=0; % only additive: 1 for multiplicative 

ops_system.ops_masses=ops_masses;

predict_horz=10;% prediction horizon
scenario_size=[2 2 1]; % branching factor

Test_points=100; % sample test points
x_rand=4*rand(ops_system.nx,Test_points)-2;

time_gpad=cell(Test_points,1);

%% Generation of tree
ops_tree.N=predict_horz; 
ops_tree.brch_ftr=ones(ops_tree.N,1);
ops_tree.brch_ftr(1:length(scenario_size))=scenario_size;
Ns=prod(ops_tree.brch_ftr);
ops_tree.nx=ops_system.nx;
ops_tree.prob=cell(ops_tree.N,1);
for i=1:ops_tree.N;
    if(i<=length(scenario_size))
        pd=rand(1,ops_tree.brch_ftr(i));
        if(i==1)
            ops_tree.prob{i,1}=pd/sum(pd);
            pm=1;
        else
            pm=pm*scenario_size(i-1);
            ops_tree.prob{i,1}=kron(ones(pm,1),pd/sum(pd));
        end
    else
        ops_tree.prob{i,1}=ones(Ns,1);
    end
end
tic
[sys_no_precond,Tree]=tree_generation_multiple(ops_system,ops_tree);
time.tree_formation=toc;
sys_no_precond.nx=ops_system.nx;
sys_no_precond.nu=ops_system.nu;
sys_no_precond.Np=ops_tree.N;
SI=scenario_index(Tree); %calculation of the scenario index.
%%
% Cost function
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
Ptree=GPAD_factor_step_smpc(sys,V,Tree);
toc
ops_GPAD.steps=200;
ops_GPAD.primal_inf=1e-3;
ops_GPAD.dual_gap=10e-3;
ops_GPAD.alpha=1/calculate_Lipschitz(sys,V,Tree);

max_size=zeros(Test_points,length(Tree.stage));

for kk=1:Test_points
    
    ops_GPAD.x0=x_rand(:,kk);
    
    % APG algorithm
    [Z_gpad_prcnd,Y_gpad_pre,time_gpad{kk}]=GPAD_algorithm_SMPC(sys,Ptree,Tree,V,ops_GPAD);
    if(~isfield(time_gpad{kk},'iterate'))
        time_gpad{kk}.iterate=ops_GPAD.steps;
    end
end

%% plots 

APG_iterations=zeros(Test_points,1);
for  i=1:Test_points
    APG_iterations(i,1)=time_gpad{i}.iterate; 
end
plot(APG_iterations)
xlabel('various test points')
ylabel('iterations')