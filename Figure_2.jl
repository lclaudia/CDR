include("DDAfunctions.jl");

SYSTEM = "ROS";
#SYSTEM = "LOR";

OUT_DIR = "CDR_LorRos"; dir_exist(OUT_DIR);

NrSyst = 2;
if SYSTEM == "ROS"
   SYS=[
       [0  0 2];  
       [0  0 3];
       [1  0 1];
       [1  0 2];
       [2  0 0];
       [2  0 3];
       [2  1 3]
       ];

   DELTA=2; dt=0.05;
   a=0.15; c=-10; b=0.2; omega1=1; omega2=1;
   MOD_par = [-omega1 -1 omega1 a b .+ randn(1,1)/100 c 1  -omega2 -1 omega2 a b .+ randn(1,1)/100 c 1];
   mm=3;  TAU=[32 9];  
end
if SYSTEM == "LOR"
   SYS=[
       [0  0 1];
       [0  0 2];
       [1  0 1];
       [1  0 2];
       [1  1 3];
       [2  0 3];
       [2  1 2]
       ];

   DELTA=1; dt=0.01;
   sigma=10; R=fill(28,NrSyst) .+ randn(NrSyst)/100; beta=8/3;
   MOD_par = repeat([-sigma sigma R[1]  -1 -1 -beta 1],NrSyst,1); MOD_par[:,3] .= R;
   MOD_par = MOD_par'[:]';
   mm=8; TAU=[6 20];
end

(MOD_nr,DIM,ODEorder,P_DDA) = make_MOD_nr(SYS,NrSyst);
TM=maximum(TAU);
X0=rand(NrSyst*DIM,1); 
CH_list1=1:3:DIM*NrSyst; 

dm=3; DDAorder=3; N_MOD=3; 
(MOD,P_DDA,SSYM) = make_MOD_new_new(N_MOD,nr_delays,DDAorder);
(MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);
PP = P_DDA[MODEL,:];

WL=3000; WS = 50; TRANS=100000;
WN=4000; LX = (WN-1)*WS + WL+TM+2*dm;

X = integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,LX,DIM,ODEorder,X0,[],CH_list1,DELTA,TRANS);

SNR_list = [-1 ; 20:-1:0];
    
NrCH=size(X,2);
CH_list = vcat(collect(combinations(1:NrCH,2))...)';
CHs  = sort(union(CH_list));
LISTE = reshape(CH_list,2,Int(length(CH_list)/2))';
LISTE = vcat(LISTE,repeat(CHs,1,2));

WN = Int(1+floor((LX-(WL+TM+2*dm))/WS));
AE=collect(0:WN-1) .* WS; AE=hcat(AE,AE .+ (WL+TM-1)); 
AE .+= 1; # julia starts with 1

NN=2;

FN_out = @sprintf("%s/CD_test_%s_%02d.jld2",OUT_DIR,SYSTEM,NN);
if !isfile(FN_out)
  CD=fill(NaN,length(SNR_list),WN,2*NN+1,size(LISTE,1),2);

  for snr=1:length(SNR_list)
    @printf("%3d\n",SNR_list[snr]);

    Y = X .* 1;
    if SNR_list[snr] >= 0
       for ch=1:2
           Y[:,ch]=add_noise(Y[:,ch],SNR_list[snr]);
       end
    end

    dY=fill(NaN,size(Y,1)-2*dm,size(Y,2));
    for ch=1:size(Y,2)
        dY[:,ch] = deriv_all(Y[:,ch],dm);
    end
    Y=Y[dm+1:end-dm,:];

    ST=fill(NaN,WN,size(Y,2),4);

    CHs=[1 2];
    M_mod=fill(NaN,WN,WL,size(PP,1),length(CHs));
    DATA =fill(NaN,WN,WL,length(CHs));
    dDATA=fill(NaN,WN,WL,length(CHs));

    for wn=1:WN
        for ch=1:length(CHs)
            data=Y[AE[wn,1]:AE[wn,2],ch]; 
            STD=std(data); data = (data .- mean(data)) ./ STD; 
            dDATA[wn,:,ch] = (dY[AE[wn,1]:AE[wn,2],ch] ./ STD)[TM+1:end];
            M_mod[wn,:,:,ch]=make_M(data,WL,TAU,PP,TM); 
            DATA[wn,:,ch] = data[TM+1:end];

            ST[wn,ch,1:3] = M_mod[wn,:,:,ch] \ dDATA[wn,:,ch]; 
            a123=M_mod[wn,:,:,ch]*ST[wn,ch,1:3];
            ST[wn,ch,4] = sqrt(mean((dDATA[wn,:,ch] .- a123).^2));
        end
    end
    ST=ST[:,:,end];
    
    for k=1:size(LISTE,1)
        ch1=LISTE[k,1];ch2=LISTE[k,2];
        for wn1=NN+1:WN-NN
            for wn2=wn1-NN:wn1+NN
                M = hcat(M_mod[wn1,:,:,ch1],M_mod[wn2,:,:,ch2]);
                A1=M*(M \ dDATA[wn1,:,ch1]); A2=M*(M \ dDATA[wn2,:,ch2]);
    
                CD[snr,wn1,wn2-wn1+(NN+1),k,1] = abs(ST[wn1,ch1] - sqrt(mean((dDATA[wn1,:,ch1] .- A1).^2)));
                CD[snr,wn1,wn2-wn1+(NN+1),k,2] = abs(ST[wn2,ch2] - sqrt(mean((dDATA[wn2,:,ch2] .- A2).^2)));
            end
        end
    end
  end
  @save FN_out CD LISTE WN
end
if isfile(FN_out)
   CD = load(FN_out,"CD");
   WN = load(FN_out,"WN");
   LISTE = load(FN_out,"LISTE");
end



NM=NN+1; # middle

#CD = load(FN_out,"CD");
cd=fill(NaN,length(SNR_list),WN,2,NrCH,NrCH);
CD[:,:,2,:,:] = dropdims(mean(CD[:,:,[NM-2:NM-1;NM+1:NM+2],:,:],dims=3),dims=3);
CD[:,:,1,:,:] = CD[:,:,NM,:,:];
CD=CD[:,:,1:2,:,:];
for l=1:size(LISTE,1)
    cd[:,:,:,LISTE[l,1],LISTE[l,2]] = CD[:,:,:,l,2];
    cd[:,:,:,LISTE[l,2],LISTE[l,1]] = CD[:,:,:,l,1];
end
CD = cd .+ 0;
cd=nothing;

idx1=diagind(ones(NrCH,NrCH));  # mit sich selbst - cyan
idx2=setdiff(1:NrCH^2,idx1);    # aehnlich        - magenta

cd = dropdims(mean(CD[:,NM+10:end-(NM-1),:,:,:],dims=2),dims=2);
cd = reshape(cd,size(cd,1),2,NrCH^2);
cd1=cd[:,2,idx1];  #  mit sich selbst
cd2=cd[:,2,idx2];  #  mit anderem
cd3=cd[:,1,idx1];  #  mit sich selbst; NM
cd4=cd[:,1,idx2];  #  mit anderem    ; NM


#scalefontsizes(2)  # only once


pp = plot(size=(800,400),margins=30*Plots.px);

plot!(pp,2:size(cd4,1),cd4[2:end,:],legend=false,c=:cyan,linewidth=2)
scatter!(pp,2:size(cd4,1),cd4[2:end,:],msw=0,legend=false,c=:cyan,markersize=4);
scatter!(pp,[-1;-1],cd4[1,:],msw=0,legend=false,c=:cyan,markersize=4);
plot!(pp,[-1;2],cd4[1:2,:],legend=false,c=:cyan,linewidth=0.5)

plot!(pp,2:size(cd2,1),cd2[2:end,:],legend=false,c=:magenta,linewidth=2)
scatter!(pp,2:size(cd2,1),cd2[2:end,:],msw=0,legend=false,c=:magenta,markersize=4);
scatter!(pp,[-1;-1],cd2[1,:],msw=0,legend=false,c=:magenta,markersize=4);
plot!(pp,[-1;2],cd2[1:2,:],legend=false,c=:magenta,linewidth=0.5)

txt=string.(replace(SNR_list, -1 => Inf));
plot!(pp,xticks=(2:5:length(SNR_list),txt[2:5:end]));
plot!(pp,xlabel="SNR",xguidefontsize=20);
plot!(pp,ylabel=L"$\mathcal{C}$",yguidefontsize=30);

plot!(pp,xlims=(-3,size(cd4,1)+1))

if SYSTEM == "ROS"
   plot!(pp,title=L"Rössler - $x$");
end
if SYSTEM == "LOR"
   plot!(pp,title=L"Lorenz - $x$");
end

display(pp)

savefig(pp,@sprintf("%s/CD_test_%s.svg",OUT_DIR,SYSTEM));
savefig(pp,@sprintf("%s/CD_test_%s.png",OUT_DIR,SYSTEM));

