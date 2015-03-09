% Where to store the header file:
%create the files which store the data
tic
ptree.d=cell(1,length(Tree.children));
for i=1:length(Tree.children)
    ptree.d{1,i}=Ptree.d{1,i}';
    if(i<=(length(Tree.children)-Ns))
        ptree.f{1,i}=Ptree.f{1,i}';
    else
        ptree.f{1,i}=Ptree.f{1,i}';
    end
end
size_of_FN = zeros(length(sys.Ft),1);
size_of_FN_cum = zeros(size(size_of_FN));
for kk=1:length(sys.Ft),
    size_of_FN(kk) = length(sys.gt{kk});
    if (kk>1), size_of_FN_cum(kk) = size_of_FN_cum(kk-1)+size_of_FN(kk-1); end
end
FN_numel = sum(size_of_FN)*sys.nx;
GN_numel = sum(size_of_FN);
P_test=Ptree.P{1};

ny=size(sys.F{1},1);
Nd=length(Tree.children);
Ns=length(Tree.leaves);
g=cell2mat(sys.g);
filedata={
    {'Data_files/GPAD_FN.h',               'Ft',           sys.Ft,                 FN_numel}
    {'Data_files/GPAD_gN.h',               'gN',          sys.gt,                 GN_numel}
    %{'Data_files/GPAD_K_GAIN.h',  'GPAD_K_GAIN',          Ptree.K,                sys.nu*sys.nx*Nd}
    %{'Data_files/GPAD_PHI.h',        'GPAD_PHI',          Ptree.Phi,              Nd*ny*sys.nu}
    %{'Data_files/GPAD_THETA.h',    'GPAD_THETA',          Ptree.Theta,            (Nd-1)*sys.nu*sys.nx+GN_numel*sys.nu}
    %{'Data_files/GPAD_D.h',            'GPAD_D',          ptree.d,                Nd*ny*sys.nx}
    %{'Data_files/GPAD_F.h',            'GPAD_F',          ptree.f,                (Nd-1)*sys.nx*sys.nx+GN_numel*sys.nx}
    {'Data_files/GPAD_Vf.h',             'V_Vf',             V.Vf,                 sys.nx*sys.nx*Ns}
    %{'Data_files/GPAD_C.h',            'GPAD_C',          Ptree.c,                 Nd*sys.nx}
    %{'Data_files/GPAD_Sigma.h',    'GPAD_SIGMA',      Ptree.sigma,                 Nd*sys.nu}
    {'Data_files/GPAD_Tree_Value.h','TREE_VALUE',       Tree.value',                 (Nd+Ns)*sys.nx}
    {'Data_files/GPAD_P.h',             'GPAD_P',            P_test,                 sys.nx*sys.nx}
    };

if(ops_system.sys_uncert)
    filedata(size(filedata,1)+1:size(filedata,1)+5,1)={
        {'Data_files/GPAD_Fc.h',                'F',            sys.F,                 ny*sys.nx*Nd}
        {'Data_files/GPAD_Gc.h',                'G',            sys.G,                 ny*sys.nu*Nd}
        {'Data_files/GPAD_g.h',                 'g',                g,                 ny*Nd}
        {'Data_files/GPAD_A.h',             'sys_A',            sys.A,                 sys.nx*sys.nx*(Nd+Ns)}
        {'Data_files/GPAD_B.h',             'sys_B',            sys.B,                 sys.nx*sys.nu*(Nd+Ns)}
        };
else
    filedata(size(filedata,1)+1:size(filedata,1)+5,1)={
        {'Data_files/GPAD_Fc.h',                'F',            sys.F{1},                 ny*sys.nx}
        {'Data_files/GPAD_Gc.h',                'G',            sys.G{1},                 ny*sys.nu}
        {'Data_files/GPAD_g.h',                 'g',            sys.g{1},                 ny}
        {'Data_files/GPAD_A.h',             'sys_A',            sys.A{1},                 sys.nx*sys.nx}
        {'Data_files/GPAD_B.h',             'sys_B',            sys.B{1},                 sys.nx*sys.nu}
        };
end
%%
for kk=1:length(filedata(:,1))
    m_matrix =filedata{kk};
    f = fopen(filedata{kk,1}{1,1},'w+');
    fprintf(f,'%d \n',filedata{kk,1}{1,4});
    m = m_matrix{3};
    if isnumeric(m), % m is a matrix
        for i=1:size(m,2),
            for j=1:size(m,1),
                fprintf(f,'%g\n', m(j,i));
            end
        end
    elseif iscell(m), % m is a cell
        for s=1:length(m),
            mm=m{s};
            for i=1:size(mm,2),
                for j=1:size(mm,1),
                    fprintf(f,'%g\n', mm(j,i));
                end
            end
        end
    end
    fclose(f);
end


