function [ S ] = system_generation(N, options )
%masses_model function is used to find the dynamic equations of N masses
%system. The no of masses is specified by N and options specify the constraints
%on the system.
% The syntax of the function is
% S=masses_generation(N,options) 
% reference is hirarchical MPC paper example.


default_options = struct('M', 0.5*ones(N,1), 'b', 0.1*0.1*ones(N,1),...
    'k',1*ones(N,1),'kp',0.1*ones(N,1), 'km', 0.1*ones(N,1), 'xmin', ...
    -5*ones(2*N,1), 'xmax', 5*ones(2*N,1), 'umin', -5*ones(N,1),'umax',...
    5*ones(N,1), 'Ts', 0.1);
ops = default_options;
if nargin==2, % User-provided options
    ops = options;
    flds = fieldnames(default_options);
    for i=1:numel(flds),
        if ~isfield(options, flds(i))
            ops.(flds{i})=default_options.(flds{i});
        end
    end
end

M=ops.M;        % M(i)=mass of body #i
b=ops.b;        % b(i)=viscous friction of body #i
k=ops.k;        % k(i)=spring of body #i
kp=ops.kp;      % kp(i)= 
km=ops.km;      % km(i)= 

% Define full A,B model
nx=2*N;
nu=N;
S.nx=nx;
S.nu=nu;
S.nr=N;


Ag=zeros(nx,nx);
Bg=zeros(nx,nu);
for i=1:nu,
    h=2*(i-1)+1;
    Ag(h,h+1)=1;  % velocity
    Ag(h+1,h+1)=-b(i)/M(i);  % friction
    Ag(h+1,h)=-k(i)/M(i);    % self-springs
    if i>1,
        Ag(h+1,h-2)=km(i)/M(i);
        Ag(h+1,h)=Ag(h+1,h)-km(i)/M(i);
    end
    if i<N,
        Ag(h+1,h+2)=kp(i)/M(i);
        Ag(h+1,h)=Ag(h+1,h)-kp(i)/M(i);
    end
    Bg(h+1,i)=1/M(i);
end

%S.A=eye(2*N)+ops.Ts*Ag;
%S.B=ops.Ts*Bg;

sysgc=ss(Ag,Bg,eye(2*N),zeros(2*N,N));

% Discretize and normalise...
sysgd=c2d(sysgc, ops.Ts);
tol=1e-6;
S.A=(1-(abs(sysgd.a)<tol)).*sysgd.a;
S.B=(1-(abs(sysgd.b)<tol)).*sysgd.b;
S.F=[eye(S.nx);-eye(S.nx);zeros(2*S.nu,S.nx)];
S.G=[zeros(2*S.nx,S.nu);eye(S.nu);-eye(S.nu)];
S.g=[ops.xmax;-ops.xmin;ops.umax;-ops.umin];

%{
nt=size(S.F,1);
for i=1:nt
    S.F(i,:)=S.F(i,:)/S.g(i);
    S.G(i,:)=S.G(i,:)/S.g(i);
    S.g(i,:)=1;
end 

S.xmin=ops.xmin;
S.xmax=ops.xmax;
S.umin=ops.umin;
S.umax=ops.umax; 
%}
end

