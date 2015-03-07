function [ L] = calculate_Lipschitz( sys,V,Tree )
%This function is used to calculate the lipschitz constant 
% of the function which is used as the step size in the gradient algorithm  
%(step 3) in the GPAD algorithm. we approximate it by an upper bound 

% consider \min_{z \in\mathcal{Z}(p)} V(z)=Z'\tilde{Q}Z+linear part 
%              such that
%              Mz<=0 [F G 0 0 
%                         F G];
% L= norm(V'\tilde{Q}V,2)<=norm(V,2)/\lambda_{\min}(\tilde{Q}};

Nd=length(Tree.stage);
Ns=length(Tree.leaves);
%M=kron(speye(Nd-Ns),[sys.F sys.G]);
if(iscell(sys.F))
    l2=norm(sys.G{1}*(V.R\sys.G{1}'),2);
    for i=2:Nd-Ns
        l2=max(l2,1/Tree.prob(i)*norm(sys.F{i}*(V.Q\sys.F{i}')+sys.G{i}*(V.R\sys.G{i}'),2));
    end
else
    l2=norm(sys.G*(V.R\sys.G'),2);
    l2=max(l2,1/min(Tree.prob(1:Nd-Ns))*norm(sys.F*(V.Q\sys.F')+sys.G*(V.R\sys.G'),2));
end 
for i=1:Ns
    l2=max(l2,1/Tree.prob(Nd-Ns+i)*norm(sys.Ft{i}*(V.Vf{i}\sys.Ft{i}')));
end 
L=l2;
end

