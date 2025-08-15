# julia


include("DDAfunctions.jl");

NrSyst=7;    
ROS=[[0  0 2];  
     [0  0 3];
     [1  0 1];
     [1  0 2];
     [2  0 0];
     [2  0 3];
     [2  1 3]
    ];
(MOD_nr,DIM,ODEorder,P) = make_MOD_nr(ROS,NrSyst);

a123= 0.21; 
a456= 0.20;
a7  = 0.18;

b1  = 0.2150;
b2  = 0.2020; 
b3  = 0.2041; 
b4  = 0.4050; 
b5  = 0.3991; 
b6  = 0.4100; 
b7  = 0.5000;

c   = 5.7;
c7  = 6.8;

MOD_par=[
         -1 -1 1 a123  b1 -c  1 
         -1 -1 1 a123  b2 -c  1
         -1 -1 1 a123  b3 -c  1
         -1 -1 1 a456  b4 -c  1
         -1 -1 1 a456  b5 -c  1
         -1 -1 1 a456  b6 -c  1
         -1 -1 1 a7    b7 -c7 1
        ];
MOD_par=reshape(MOD_par',size(ROS,1)*NrSyst)';

FromTo=[ 
        [3 0  0 1   7 0  0 1];                                        # coupling
        [7 0  0 1   1 0  0 1];
        [7 0  0 1   4 0  0 1];
        [7 0  0 1   5 0  0 1]
        ];
PF = "__FirstExample";                                                # parts of file names

#=                                                                    # unconnected network
FromTo=[];
PF = "__Empty";
=#

II=make_MOD_nr_Coupling(FromTo,DIM,P);                                # MOD_nr part for coupling
epsilon=0.15;                                                         # coupling strength
MOD_par_add=repeat([epsilon -epsilon],size(FromTo,1),1)'[:]';         # MOD_par for coupling part

TAU=[32 9]; TM=maximum(TAU); dm=4;                                    # DDA parameters
WL=4000;WS=2000;                                                      # window length and window shift for DDA 
WN=100;                                                               # assign window number 
LL=WS*(WN-1)+WL+TM+2*dm;                                              # ajust integration length

TRANS=20000;                                                          # transient
dt=0.05;                                                              # integration step size
X0=rand(DIM*NrSyst,1);                                                # initial conditions
DATA_DIR="DATA"; dir_exist(DATA_DIR);                                 # DATA folder
noise="NoNoise";NOISE="NoNoise";                                      # parts of file names
FN=@sprintf("%s%sCD_DDA_data_%s__WL%d_WS%d_WN%d%s.ascii", 
            DATA_DIR,SL,noise,WL,WS,WN,PF);                           # noise free data file
CH_list=1:DIM:DIM*NrSyst;                                             # only x
DELTA=2;                                                              # every second data point
if !isfile(FN)
   integrate_ODE_general_BIG([MOD_nr II],[MOD_par MOD_par_add],       # encoding of the coupled systems
                             dt,                                      # step size of num. integration
                             LL,                                      # length 
                             DIM*NrSyst,ODEorder,X0,                  # parameters
                             FN,
                             CH_list,DELTA,
                             TRANS);
end

SNRadd_list= 20:-1:0;

MakeDataNoise(PF,noise,SNRadd_list);                                  # add noise


DDA_DIR="DDA"; dir_exist(DDA_DIR);                                    # DDA folder

nr_delays=2; 
DDAmodel=[[0 0 1];  
          [0 0 2]; 
          [1 1 1]];
(MODEL, L_AF, DDAorder)=make_MODEL(DDAmodel);                         # DDA model

NrCH=NrSyst; CH=collect(1:NrCH);  
LIST=collect(combinations(CH,2)); 
LL1=vcat(LIST...)';
LIST=reduce(hcat,LIST)';                                              # pairwise combinations

RunDDA(PF,NOISE,SNRadd_list);                                         # run DDA


(C,E)=makeCE(PF,NOISE,SNRadd_list);                                   # read outputs

c=reshape(C[:,:,1],NrSyst^2); IDX=reverse(sortperm(c[:])); IDX=IDX[1:NrSyst^2-NrSyst];
c=reshape(C,NrSyst^2,length(SNRadd_list)+1); c=c[IDX,:];
e=reshape(E,NrSyst^2,length(SNRadd_list)+1); e=e[IDX,:];

S=[repeat(1:NrSyst,NrSyst,1) repeat(1:NrSyst,1,NrSyst)'[:]][IDX,:];

p1=Plots.palette(cgrad(:cool,scale=:log,rev=true),42);
p2=Plots.palette(cgrad(:cool,scale=:log),42);

A=reshape(C[:,:,1],NrSyst^2); idx=reverse(sortperm(A[:])); A=A[idx[1:NrSyst^2-NrSyst]]; 
A=A .- A[end]; A=A ./ A[1] .* 1000; A=Int.(floor.(A)); A[A .== 0] .= 1;
p1=Plots.palette(cgrad(:cool,scale=:log),1000); p1=p1[A];

pp=plot(size=(1000,550),margins=10*Plots.px,layout=@layout[a{0.8w} [b;c]]);

scatter!(pp,subplot=1,c[:,1:end]',e[:,1:end]',msw=0,palette=p1,markersize=5,xscale=:log10,yscale=:log10,label="")
plot!(pp,subplot=1,c',e',msw=0,palette=p1,linewidth=0.2,xscale=:log10,yscale=:log10,label="",grid=false)

plot!(pp,subplot=1,xlabel=L"$\mathcal{C}$");
plot!(pp,subplot=1,ylabel=L"$\mathcal{E}$");

#=
for k=1:size(FromTo,1)+3
    txt=join([string(S[k,1]) " \u21FE " string(S[k,2])])
    annotate!(pp,subplot=1,c[k,1],e[k,1],Plots.text(txt,18,p1[k]));
end
=#
display(pp)

cc=C[:,:,1]; cc=cc .- minimum(cc); cc=cc ./ maximum(cc);
cc[diagind(cc)] .= NaN;

heatmap!(pp,subplot=2,cc,aspect_ratio=:equal,c=p2,colorbar=false,grid=false,axis=([], false))
display(pp)


GR.setarrowsize(0.5);

MS = [1,1,1,2,2,2,3];
colors = [colorant"plum2", colorant"mistyrose1", colorant"lavender"];

A=C[:,:,1];
A[A .== 0] .= 1;
A = A .- minimum(A);
A[A .== maximum(A)] .= 0;
A = A ./ maximum(A);

A[A .< 0.33] .= 0;

graphplot!(pp,subplot=3,A,
              method=:circular,nodeshape=:circle,
              names=1:7,
              markersize=0.15,
              fontsize=20,
              linewidth=3,
              linealpha=1,
              markercolor = colors[MS],
              nodestrokecolor=colors[MS],
              arrow=arrow(:closed,10),
              )

display(pp)

############

pp1=plot(size=(1000,200),margins=10*Plots.px,layout=(1,5));
k=1;cc=C[:,:,1]; M1=minimum(cc); cc=cc .- M1; M2=maximum(cc); cc=cc ./ M2;
cc[diagind(cc)] .= NaN;
heatmap!(pp1,subplot=Int((k - 1)/5 + 1),cc,
             aspect_ratio=:equal,c=p2,colorbar=false,
             grid=false,axis=([], false),clims=(0,1),
             title="no noise")

for k=6:5:22
    cc=C[:,:,k]; cc=cc .- M1; cc=cc ./ M2; cc[diagind(cc)] .= NaN;
    heatmap!(pp1,subplot=Int((k - 1)/5 + 1),cc,
             aspect_ratio=:equal,c=p2,colorbar=false,
             grid=false,axis=([], false),clims=(0,1),
             title=@sprintf("SNR = %d dB",SNRadd_list[k]))
end
display(pp1)
