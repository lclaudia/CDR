### run Figure_3.jl first


TEn=fill(NaN,WN,22,2);
for SNR = 0:20
    FN = @sprintf("../EXP_RosPalus07/TE__SNR_%02d.dat",SNR);
    TEn[:,SNR+1,:] = readdlm(FN);
end
TEn[:,22,:] = TE;

GC=fill(NaN,WN,22,2);
for SNR = 0:20
    FN = @sprintf("../EXP_RosPalus07/GC__SNR_%02d.dat",SNR);
    GC[:,SNR+1,:] = readdlm(FN);
end
GC[:,22,:] = Granger;

CDnoise = [CDnoise CD];




f_list=[-0.15,0.05,0.135,0.2];

p1=Plots.palette(cgrad(:jet,rev=true),8);
txt=string.(20:-5:0);

pp=plot(size=(800,1800),layout=@layout[a;b;c],left_margin=50*Plots.px, bottom_margin=25*Plots.px)

for k=1:length(f_list)
    f=findfirst(epsilon .>= f_list[k]);

    cd=dropdims(mean(CDnoise[f .+ (-20:20),:,:],dims=1),dims=1); 
    gc=dropdims(mean(GC[f .+ (-20:20),:,:],dims=1),dims=1);
    te=dropdims(mean(TEn[f .+ (-20:20),:,:],dims=1),dims=1);

    plot!(pp,subplot=1,log10.(cd[end:-1:1,[2;1]]),palette=p1,label="",linewidth=1);
    plot!(pp,subplot=2,log10.(gc[end:-1:1,:]),msw=0,palette=p1,label="",linewidth=1);
    plot!(pp,subplot=3,te[end:-1:1,:],palette=p1,label="",linewidth=1);

    display(pp)
end


for k=1:length(f_list)
    f=findfirst(epsilon .>= f_list[k]);

    cd=dropdims(mean(CDnoise[f .+ (-20:20),:,:],dims=1),dims=1); 
    gc=dropdims(mean(GC[f .+ (-20:20),:,:],dims=1),dims=1);
    te=dropdims(mean(TEn[f .+ (-20:20),:,:],dims=1),dims=1);

    scatter!(pp,subplot=1,log10.(cd[end:-1:1,[2;1]]),msw=0,palette=p1,label="");
    scatter!(pp,subplot=2,log10.(gc[end:-1:1,:]),msw=0,palette=p1,label="");
    scatter!(pp,subplot=3,te[end:-1:1,:],msw=0,palette=p1,label="");

    display(pp)
end



plot!(pp,subplot=1,ylabel=L"$\mathcal{C}$",yguidefontsize=20);
plot!(pp,subplot=2,ylabel="GC",yguidefontsize=20);
plot!(pp,subplot=3,ylabel="TE",yguidefontsize=20);

for k=1:3
    plot!(pp,subplot=k,xticks=(2:5:22,txt));
    plot!(pp,subplot=k,xlabel="SNR",xguidefontsize=20);
end

display(pp)



savefig(pp,"Palus_CDR.svg");
