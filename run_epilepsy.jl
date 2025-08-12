# julia -p10
# /cnl/xaos/julia-1.11.1/bin/julia -p50

@everywhere MacSalk="SALK";

if MacSalk=="SALK"
   @everywhere pfad_home="/home/claudia/";
end
if MacSalk=="MAC"
   @everywhere pfad_home="/Users/claudia/";
end

@everywhere include(join([pfad_home "/TOOLS/DDAfunctions.jl"]));

@everywhere SL="/";

EDF_file = "patient1_S05__01_03.edf";

NrCH=78;
CH=1:NrCH;
ONSET_CH = sort([45:52; 85:92; 61:68; 101:108] .- 30);

TAU=[7 10]; TM = maximum(TAU); dm=4;
WL=500; WS=50;

nr_delays=2; 
DDAmodel=[[0 0 0 1];  
          [0 0 0 2]; 
          [1 1 1 1]];
(MODEL, L_AF, DDAorder)=make_MODEL(DDAmodel);                         # DDA model

DDA_DIR="DDA_patient1"; dir_exist(DDA_DIR);

LIST=reduce(hcat,collect(combinations(CH,2)))';                       # original data
LIST=vcat(LIST, LIST .+ NrCH);                                        # data with added noise
DELTA=20;
N=Int64(ceil(size(LIST,1)/DELTA));

FN_DATA = EDF_file;
FN_DDA  = @sprintf("%s%s%s",DDA_DIR,SL,replace(EDF_file,".edf" => ".DDA"));
FN_ALL = @sprintf("%s%s%s",DDA_DIR,SL,replace(EDF_file,".edf" => ".jld2"));
if !isfile(FN_ALL)
   @printf("%s\n",FN_ALL);

   if !isfile(join([FN_DDA,"_ST"]))
      if Sys.iswindows()
         if !isfile("run_DDA_ASCII.exe")
            run(`cp run_DDA_ASCII run_DDA_ASCII.exe`);
         end

         CMD=".\\run_DDA_ASCII.exe";
      else
         CMD="./run_DDA_ASCII";
      end
   
      CMD = "$CMD -EDF";                                 
      CMD = "$CMD -MODEL $(join(MODEL," "))"                    
      CMD = "$CMD -TAU $(join(TAU," "))"                        
      CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays"    
      CMD = "$CMD -DATA_FN $FN_DATA -OUT_FN $FN_DDA"        
      CMD = "$CMD -WL $WL -WS $WS";                            
      CMD = "$CMD -SELECT 1 0 0 0";                               # ST-DDA               
    
      if Sys.iswindows()
         run(Cmd(string.(split(CMD, " "))));
      else
         run(`sh -c $CMD`);
      end

      rm(@sprintf("%s.info",FN_DDA));     
   end
   
   @sync @distributed for n_N=1:N
       FN_DDAn=@sprintf("%s%s%s__%03d.DDA",DDA_DIR,SL,replace(EDF_file,".edf" => ""),n_N);
       
       n=collect(1:DELTA) .+ (n_N-1)*DELTA; n=n[n.<=size(LIST,1)];
       LL1=LIST[n,:]; LL1=vcat(LL1'...)';
   
       if !isfile(join([FN_DDAn,"_CT"]))
          if Sys.iswindows()
             if !isfile("run_DDA_ASCII.exe")
                run(`cp run_DDA_ASCII run_DDA_ASCII.exe`);
             end
 
             CMD=".\\run_DDA_ASCII.exe";
          else
             CMD="./run_DDA_ASCII";
          end

          CMD = "$CMD -EDF";                                 
          CMD = "$CMD -MODEL $(join(MODEL," "))"                    
          CMD = "$CMD -TAU $(join(TAU," "))"                        
          CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays"    
          CMD = "$CMD -DATA_FN $FN_DATA -OUT_FN $FN_DDAn"        
          CMD = "$CMD -WL $WL -WS $WS";                            
          CMD = "$CMD -SELECT 0 1 0 0"                            # CT-DDA                           
          CMD = "$CMD -CT_CH_list $(join(LL1," "))";              # all pairwise channels
          CMD = "$CMD -WL_CT 2 -WS_CT 2";
            
          if Sys.iswindows()
             run(Cmd(string.(split(CMD, " "))));
          else
             run(`sh -c $CMD`);
          end
  
          rm(@sprintf("%s.info",FN_DDAn));     
       end
   end
   
   @sync @distributed for n_N=1:N
      FN_DDAn=@sprintf("%s%s%s__%03d.DDA",DDA_DIR,SL,replace(EDF_file,".edf" => ""),n_N);
       
       n=collect(1:DELTA) .+ (n_N-1)*DELTA; n=n[n.<=size(LIST,1)];
       LL1=LIST[n,:]; LL1=vcat(LL1'...)';

       if !isfile(join([FN_DDAn,"_CD_DDA_ST"]))
          if Sys.iswindows()
             if !isfile("run_DDA_ASCII.exe")
                run(`cp run_DDA_ASCII run_DDA_ASCII.exe`);
             end
 
             CMD=".\\run_DDA_ASCII.exe";
          else
             CMD="./run_DDA_ASCII";
          end
   
          CMD = "$CMD -EDF";                                 
          CMD = "$CMD -MODEL $(join(MODEL," "))"                    
          CMD = "$CMD -TAU $(join(TAU," "))"                        
          CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays"    
          CMD = "$CMD -DATA_FN $FN_DATA -OUT_FN $FN_DDAn"        
          CMD = "$CMD -WL $WL -WS $WS";                            
          CMD = "$CMD -SELECT 0 0 1 0";                           # CD-DDA
          CMD = "$CMD -PAIRS $(join(LL1," "))";                   # all pairwise channels
        
          if Sys.iswindows()
             run(Cmd(string.(split(CMD, " "))));
          else
             run(`sh -c $CMD`);
          end
  
          rm(@sprintf("%s.info",FN_DDAn));     
       end
   end

   ST=readdlm(join([FN_DDA,"_ST"])); 
   T=ST[:,1:2]; ST=ST[:,3:end];
   WN=size(T,1);

   #=   plot a1
   a1=ST[:,1:4:NrCH*4];
   heatmap(a1[:,setdiff(1:NrCH,54)]',legend=false,c=:jet,clims=(-0.25,0.4))
   =#
   
   rhoS=ST[:,L_AF:L_AF:end];
   ST = nothing; GC.gc();
   
   E=fill(NaN,WN,NrCH*2,NrCH*2);
   for n_N=1:N
       @printf("%3d ",n_N)
       FN_DDAn=@sprintf("%s%s%s__%03d.DDA",DDA_DIR,SL,replace(EDF_file,".edf" => ""),n_N);
   
       n=collect(1:DELTA) .+ (n_N-1)*DELTA; n=n[n.<=size(LIST,1)];
       LL1=LIST[n,:] .- CH[1] .+ 1;
   
       CT=readdlm(join([FN_DDAn,"_CT"]));
       CT=CT[:,3:end];
       CT=CT[:,L_AF:L_AF:end];
   
       for l=1:size(LL1,1)
           ch1=LL1[l,1];ch2=LL1[l,2];
           E[:,ch1,ch2] = abs.( dropdims(mean(rhoS[:,[ch1,ch2]],dims=2),dims=2) ./ CT[:,l] .- 1 );
           E[:,ch2,ch1] = E[:,ch1,ch2];
       end
       CT = nothing; GC.gc();
   end
   @printf("\n");
       
   C=fill(NaN,WN,NrCH*2,NrCH*2);
   for n_N=1:N
       @printf("%3d ",n_N)
       FN_DDAn=@sprintf("%s%s%s__%03d.DDA",DDA_DIR,SL,replace(EDF_file,".edf" => ""),n_N);
   
       n=collect(1:DELTA) .+ (n_N-1)*DELTA; n=n[n.<=size(LIST,1)];
       LL1=LIST[n,:] .- CH[1] .+ 1;
   
       CD=readdlm(join([FN_DDAn,"_CD_DDA_ST"]));
       CD=CD[:,3:end];
       CD=reshape(CD,WN,2,size(LL1,1));
   
       for l=1:size(LL1,1)
           ch1=LL1[l,1];ch2=LL1[l,2];
           C[:,ch1,ch2] = CD[:,2,l];
           C[:,ch2,ch1] = CD[:,1,l];
       end
       CD = nothing; GC.gc();
   end
   @printf("\n\n");

   @save FN_ALL C E rhoS T WN
   E = nothing; C = nothing; GC.gc();
end

@load FN_ALL C E T WN ;
C1 = C[:,1:NrCH,1:NrCH]; 
E1 = E[:,1:NrCH,1:NrCH]; 
C2 = C[:,NrCH+1:end,NrCH+1:end]; 
E2 = E[:,NrCH+1:end,NrCH+1:end]; 
E = nothing; C = nothing; GC.gc();





eLABEL=["LFP"; "LCG"; "LAT"; "LMT"; "LHP"; "LOC"; "LTH"; "LSU"; 
        "RFP"; "RCG"; "RAT"; "RMT"; "RHP"; "ROC"; "RTH"; "RSU"]

e_list = [
       [34:38; 15:22; 25:30],      #LFP

       [31:33; 23; 24],            #LCG

       [3:5; 10:14],               #LAT
       Int[],  
       [1:2; 7:9],                 #LHP
       Int[],  
       Int[],  
       Int[], 
          
       [55:62; 71:78; 63; 65:70],  #RFP
       Int[],  
       [40:44; 50:52],             #RAT
       Int[],  
       [39; 47:49],                #RHP
       Int[],  
       Int[],
       [64]                        #RSU
       ];

e_NotZero = findall(x -> x == 1, length.(e_list)' .!= 0 );
e_NotZero = [i[2] for i in e_NotZero];

SR=500;
t=(T[:,1] .+ 1 .+ TM .+ dm) ./ SR ./ 60;

SEQ=1:length(e_NotZero);

CHs=vcat(e_list[e_NotZero][SEQ][:]...);
L_e_list=map(x -> length(e_list[e_NotZero][x]),SEQ);

IND=setdiff(1:length(CHs)^2,diagind(C1[1,CHs,CHs]));
c1=reshape(C1[:,CHs,CHs],WN,length(CHs)^2)[:,IND];
c2=reshape(C2[:,CHs,CHs],WN,length(CHs)^2)[:,IND];
e1=reshape(E1[:,CHs,CHs],WN,length(CHs)^2)[:,IND];
e2=reshape(E2[:,CHs,CHs],WN,length(CHs)^2)[:,IND];

c2[c2 .< 0.01] .= NaN;
c1[c1 .< 0.01] .= NaN;

BETA=asin.((c2 .- c1) ./ sqrt.( (c1 .- c2).^2 + (e1 .- e2).^2) );

SG = plot(size=(1000,1000),layout=(2,2));

heatmap!(SG,subplot=1,
         t,1:length(IND),c1',
         c=:jet,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,clims=(0.1,0.06)
         )
hline!(SG,subplot=1,[cumsum(L_e_list[:]) .* (length(CHs)-1) .+ 0.5],legend=false,c=:black,linewidth=2)
Y = vcat(0,cumsum(L_e_list[:]) .* (length(CHs)-1));
Y = (Y .- [0;diff(Y)./2])[2:end];
heatmap!(SG,subplot=1,yticks=(Y,eLABEL[e_NotZero][SEQ]),xlabel="time [min]",clims=(0.01,0.06))
display(SG)


heatmap!(SG,subplot=2,
         t,1:length(IND),c2',
         c=:jet,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,clims=(0.1,0.06)
         )
hline!(SG,subplot=2,[cumsum(L_e_list[:]) .* (length(CHs)-1) .+ 0.5],legend=false,c=:black,linewidth=2)
heatmap!(SG,subplot=2,yticks=(Y,eLABEL[e_NotZero][SEQ]),xlabel="time [min]",clims=(0.01,0.06))
display(SG)


ALPHA = BETA .* 1;
ALPHA[ALPHA .>= 0] .= 0;
ALPHA[ALPHA .< 0]  .= 1;
ALPHA[ALPHA .== 0] .= NaN;

heatmap!(SG,subplot=3,
         t,1:length(IND),(c1 .* ALPHA)',
         c=:jet,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,clims=(0.1,0.06)
         )
hline!(SG,subplot=3,[cumsum(L_e_list[:]) .* (length(CHs)-1) .+ 0.5],legend=false,c=:black,linewidth=2)
heatmap!(SG,subplot=3,yticks=(Y,eLABEL[e_NotZero][SEQ]),xlabel="time [min]",clims=(0.01,0.06))
display(SG)



ALPHA = BETA .* 1;
ALPHA[ALPHA .<= 0] .= 0;
ALPHA[ALPHA .> 0]  .= 1;
ALPHA[ALPHA .== 0] .= NaN;

heatmap!(SG,subplot=4,
         t,1:length(IND),(c1 .* ALPHA)',
         c=:jet,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,clims=(0.1,0.06)
         )
hline!(SG,subplot=4,[cumsum(L_e_list[:]) .* (length(CHs)-1) .+ 0.5],legend=false,c=:black,linewidth=2)
heatmap!(SG,subplot=4,yticks=(Y,eLABEL[e_NotZero][SEQ]),xlabel="time [min]",clims=(0.01,0.06))
display(SG)


