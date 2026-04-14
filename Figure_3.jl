include("DDAfunctions.jl");

SYSTEM="RosPalus07";

Sdir=@sprintf("EXP_%s",SYSTEM);
#=
if isdir(Sdir)
   rm(Sdir, recursive=true)
end;
=#
DATA_DIR=Sdir; dir_exist(DATA_DIR);

NrSyst = 2;
   ROS=[
       [0  0 2];  
       [0  0 3];
       [1  0 1];
       [1  0 2];
       [2  0 0];
       [2  0 3];
       [2  1 3]
       ];
(MOD_nr,DIM,ODEorder,P) = make_MOD_nr(ROS,NrSyst);

TRANS=100000;
DELTA=2; dt=0.05;
a=0.15; c=-10; b=0.2; omega1=1.015; omega2=0.985;

MOD_par = [-omega1 -1 omega1 a b c 1  -omega2 -1 omega2 a b c 1];

FromTo=[ 
        [ 1 3  0 1 ];  
        [ 1 3  0 4 ];  
        ];
II = make_CouplingMOD_nr(FromTo,DIM,P);
# P[II[2:2:end] .+ 1,:]

MOD_par_add = [1 -1];

AddTerms=length(MOD_par_add);

FN_data = @sprintf("%s%sdata_chirp.ascii",DATA_DIR,SL); 
start = -0.25; delta=0.0000001; ende=0.25;
X0=rand(NrSyst*DIM,1); CH_list1=1:3:DIM*NrSyst; 

mm=3;  TAU=[32 9]; 
TM=maximum(TAU);

if !isfile(FN_data)
  CMD="./i_ODE_epsilon_chirp";

  BUF = "-MODEL $(join([MOD_nr II], " "))";
  BUF = "$BUF -PAR $(join([MOD_par MOD_par_add], " "))";
  BUF = "$BUF -ANF $(join(X0," "))";
  BUF = "$BUF -dt $dt";
  BUF = "$BUF -DIM $DIM";
  BUF = "$BUF -order $ODEorder";
  BUF = "$BUF -TRANS $TRANS";
  BUF = "$BUF -FILE $FN_data";
  BUF = "$BUF -CH_list $(join(CH_list1," ")) -DELTA $DELTA";
  BUF = "$BUF -epsilon $(@sprintf("%d %20.15f %20.15f %20.15f",AddTerms,start,delta,ende))";

  run(pipeline(IOBuffer(BUF),`xargs $CMD`));
end

dm=3; DDAorder=3; N_MOD=3; 
(MOD,P_DDA,SSYM) = make_MOD_new_new(N_MOD,nr_delays,DDAorder);

(MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);
PP = P_DDA[MODEL,:];

DDA_DIR = Sdir; dir_exist(DDA_DIR);

X = readdlm(FN_data); 

WL = 3000; WS=100; 
XY_list=nothing; SELECT=[0 0 1 0 1 1 1]; CH_list=[1 2];

#=
#### run_TE_Palus.m  ---  takes a long time
#### run_GC_Palus.m
=#

Granger=readdlm("../EXP_RosPalus07/GC_data_chirp.dat"); 
Granger[Granger .<= 0] .=0;

TE=readdlm("../EXP_RosPalus07/TE_data_chirp.dat"); 

FN_noise_pairs = @sprintf("%s%snoise_pairs.jld2",DATA_DIR,SL);
if !isfile(FN_noise_pairs)
   CH_list=[1 2];
   (ST, CT, CD, E, SY, RC, WN, AE)=make_SELECT(X,dm,WL,WS,TAU,PP,TM,TAU,PP,CH_list,XY_list,SELECT);
   (THETA, WN, AE) = make_THETA(X,dm,WL,WS,TM,CH_list);

   Y = fill(NaN,size(X,1),2,length(0:20));
   for SNR = 0:20
       for k=1:2
           Y[:,k,SNR+1] = add_noise(X[:,k],SNR);
       end
   end
   for SNR = 0:20
       FN_data = @sprintf("%s%sdata_chirp__SNR_%02d.ascii",DATA_DIR,SL,SNR); 
       fid=open(FN_data,"w");
       for k=1:size(X,1)
           @printf(fid,"%15.10lf %15.10lf\n",Y[k,1,SNR+1],Y[k,2,SNR+1]);
       end
       close(fid);
   end
       
   Y=reshape(Y,size(X,1),2*length(0:20));
   
   CH_list = collect(1:2*length(0:20))'
   (ST, CT, CDnoise, Enoise, SY, RC, WN, AE)=make_SELECT(Y,dm,WL,WS,TAU,PP,TM,TAU,PP,CH_list,XY_list,SELECT);
   (THETAnoise, WN, AE) = make_THETA(Y,dm,WL,WS,TM,CH_list);
   
   @save FN_noise_pairs CDnoise Enoise THETAnoise CD E THETA
end
if isfile(FN_noise_pairs)
   @load "EXP_RosPalus07/noise_pairs.jld2"
end

WN = Int(1+floor((size(X,1)-(WL+TM+2*dm))/WS));
epsilon=range(start,ende,WN);
ALPHA=0.05;

TXT=["D \u21FE R";"R \u21FE D"];
TXT=["";""];

#scalefontsizes(2)

pp=plot(size=(1000,1400),layout=(5,1),left_margin=50*Plots.px, bottom_margin=25*Plots.px)

scatter!(pp,subplot=1,epsilon,CD[:,:,2],msw=0,label=TXT[1],markersize=2,c=:blue,alpha=ALPHA)
scatter!(pp,subplot=1,epsilon,CD[:,:,1],msw=0,label=TXT[2],markersize=2,c=:red,alpha=ALPHA)
scatter!(pp,subplot=2,epsilon,THETA,msw=0,label="",markersize=2,c=:black,alpha=ALPHA)
scatter!(pp,subplot=3,epsilon,E,msw=0,label="",markersize=2,c=:black,alpha=ALPHA)

scatter!(pp,subplot=4,epsilon,log10.(Granger[:,1]),msw=0,label=TXT[1],markersize=2,c=:blue,alpha=ALPHA)
scatter!(pp,subplot=4,epsilon,log10.(Granger[:,2]),msw=0,label=TXT[2],markersize=2,c=:red,alpha=ALPHA)
scatter!(pp,subplot=4,ylims=(-5,1));

scatter!(pp,subplot=5,epsilon,TE[:,1],msw=0,label=TXT[1],markersize=2,c=:blue,alpha=ALPHA)
scatter!(pp,subplot=5,epsilon,TE[:,2],msw=0,label=TXT[2],markersize=2,c=:red,alpha=ALPHA)

plot!(pp,subplot=1,ylabel=L"$\mathcal{C}$",yguidefontsize=30);
plot!(pp,subplot=2,ylabel=L"$\mathcal{O}$",yguidefontsize=30);
plot!(pp,subplot=3,ylabel=L"$\mathcal{E}$",yguidefontsize=30);
plot!(pp,subplot=4,ylabel="GC",yguidefontsize=30);
plot!(pp,subplot=5,ylabel="TE",yguidefontsize=30);

plot!(pp,subplot=5,xlabel=L"$\epsilon$",xguidefontsize=30);

for k=1:5
    plot!(pp,subplot=k,xlims=(start,ende));
    vline!(pp,subplot=k,[-0.041, 0.12, 0, 0.1, 0.154],label="",c=:black,linewidth=2)
    vline!(pp,subplot=k,[-0.15,0.05,0.135,0.2],label="",c=:magenta,linewidth=2)
end

display(pp)

savefig(pp,"Palus_COE.png");

