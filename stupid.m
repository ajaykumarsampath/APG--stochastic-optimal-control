%% 
% Copyright 2014, Gurobi Optimization, Inc.

% This example formulates and solves the following simple QP model:
%  minimize
%      x^2 + x*y + y^2 + y*z + z^2 + 2 x
%  subject to
%      x + 2 y + 3 z >= 4
%      x +   y       >= 1
%
% It solves it once as a continuous model, and once as an integer
% model.

clear model;

Np=1;
nc=size(sys_no_precond.F,1);
nc_t=size(sys_no_precond.Ft{1},1);
krow=nc*Np+nc_t;
nvar=sys_no_precond.nx+sys_no_precond.nu;
model.Q=sparse(Np*nvar+sys_no_precond.nx,Np*nvar+sys_no_precond.nx);
model.A=sparse(krow+(Np+1)*sys_no_precond.nx,nvar*Np+sys_no_precond.nx);
for i=1:Np
    for j=1:sys_no_precond.nx
        model.varnames{(i-1)*nvar+j}=['x_' num2str(i),num2str(j)];
    end
    for j=1:sys_no_precond.nu
        model.varnames{(i-1)*nvar+sys_no_precond.nx+j}=['u_' num2str(i),num2str(j)];
    end
    model.Q((i-1)*nvar+1:i*nvar,(i-1)*nvar+1:i*nvar) = speye(sys_no_precond.nx+sys_no_precond.nu);
    model.A((i-1)*nc+1:i*nc,(i-1)*nvar+1:i*nvar)=sparse([sys_no_precond.F sys_no_precond.G]);
    model.A(krow+(i-1)*sys_no_precond.nx+1:krow+i*sys_no_precond.nx,(i-1)*nvar+1:i*nvar+sys_no_precond.nx)...
        =sparse([sys_no_precond.A{i} sys_no_precond.B{i} -speye(sys_no_precond.nx)]);
end
for j=1:sys_no_precond.nx
    model.varnames{Np*nvar+j}=['x_' num2str(Np+1),num2str(j)];
end

model.A(Np*nc+1:Np*nc+nc_t,Np*nvar+1:Np*nvar+sys_no_precond.nx)=sparse(sys_no_precond.Ft{1});
model.A(krow+Np*sys_no_precond.nx+1:krow+(Np+1)*sys_no_precond.nx,1:sys_no_precond.nx)=speye(sys_no_precond.nx);

model.rhs=zeros(krow+(Np+1)*sys_no_precond.nx,1);
model.rhs(1:nc*Np,1)=kron(ones(Np,1),sys_no_precond.g);
model.rhs(nc*Np+1:krow,1)=sys_no_precond.gt{1}';
model.sense=[repmat('<',krow,1);repmat('=',(Np+1)*sys_no_precond.nx,1)];
model.rhs(end-sys_no_precond.nx+1:end)=rand(sys_no_precond.nx,1);
model.obj = zeros(1,nvar*Np+sys_no_precond.nx);
model.lb=-100*ones(Np*nvar+sys_no_precond.nx,1);
gurobi_write(model, 'qp.lp');

results = gurobi(model);

%%

clear model;

Np=1;
nc=size(sys_no_precond.F,1);
nc_t=size(sys_no_precond.Ft{1},1);
krow=nc*Np+nc_t;
nvar=sys_no_precond.nx+sys_no_precond.nu;
model.Q=sparse(Np*nvar+sys_no_precond.nx,Np*nvar+sys_no_precond.nx);
model.A=sparse((Np+1)*sys_no_precond.nx,nvar*Np+sys_no_precond.nx);
for i=1:Np
    for j=1:sys_no_precond.nx
        model.varnames{(i-1)*nvar+j}=['x_' num2str(i),num2str(j)];
    end
    for j=1:sys_no_precond.nu
        model.varnames{(i-1)*nvar+sys_no_precond.nx+j}=['u_' num2str(i),num2str(j)];
    end
    model.Q((i-1)*nvar+1:i*nvar,(i-1)*nvar+1:i*nvar) = speye(sys_no_precond.nx+sys_no_precond.nu);
    model.A((i-1)*sys_no_precond.nx+1:i*sys_no_precond.nx,(i-1)*nvar+1:i*nvar+sys_no_precond.nx)...
        =sparse([sys_no_precond.A{i} sys_no_precond.B{i} -speye(sys_no_precond.nx)]);
end
for j=1:sys_no_precond.nx
    model.varnames{Np*nvar+j}=['x_' num2str(Np+1),num2str(j)];
end

model.A(Np*sys_no_precond.nx+1:(Np+1)*sys_no_precond.nx,1:sys_no_precond.nx)=speye(sys_no_precond.nx);

model.rhs=zeros((Np+1)*sys_no_precond.nx,1);
model.sense=repmat('=',(Np+1)*sys_no_precond.nx,1);
model.rhs(end-sys_no_precond.nx+1:end)=rand(sys_no_precond.nx,1);
model.obj = zeros(1,nvar*Np+sys_no_precond.nx);
model.lb=
gurobi_write(model, 'qp1.lp');

results = gurobi(model);

%%
model.Q = speye(sys_no_precond.nx+sys_no_precond.nu);
model.A =sparse(2*sys_no_precond.nx+sys_no_precond.nu); 
model.A(1:nc,1:sys_no_precond.nx+sys_no_precond.nu)=sparse([sys_no_precond.F sys_no_precond.G]);
model.A(nc+1:nc+sys_no_precond.nx,1:sys_no_precond.nx)=speye(sys_no_precond.nx);
model.obj = zeros(1,sys_no_precond.nx+sys_no_precond.nu);
model.rhs(1:nc) =sys_no_precond.g;
model.rhs(nc+1:nc+sys_no_precond.nx)=rand(sys_no_precond.nx,1);
model.sense = [repmat('<',nc,1);repmat('=',sys_no_precond.nx,1)];



for v=1:length(names)
    fprintf('%s %e\n', names{v}, results.x(v));
end

fprintf('Obj: %e\n', results.objval);

model.vtype = 'B';

results  = gurobi(model);

for v=1:length(names)
    fprintf('%s %e\n', names{v}, results.x(v));
end

fprintf('Obj: %e\n', results.objval);
