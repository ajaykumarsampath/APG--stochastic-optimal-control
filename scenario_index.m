function SI = scenario_index(Tree)

% SCENARIO_INDEX compute scenario indices such that SI{j} is an array with
% all node indices related to scenario j.
%
% by D. Bernardini

ns = numel(Tree.leaves); % no. of scenarios
SI = cell(numel(Tree.leaves),1); % initialization
for s=1:ns
    SI{s} = zeros(Tree.stage(Tree.leaves(s))+1,1);
    SI{s}(end) = Tree.leaves(s);
    for k=Tree.stage(Tree.leaves(s)):-1:1
        SI{s}(k) = Tree.ancestor(SI{s}(k+1));
    end
end