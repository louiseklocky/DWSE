function output=nSE(popu,fit,input)
FES=input.FES;
popuSize=input.popusize;
D=input.D;
lu=input.lu;
f=input.f;
bsf_fit_var=input.bsf_fit_var;
optimalChart=input.optimalChart;
benchmark=input.benchmark;
optimalChart_interval=input.optimalChart_interval;
run_id=input.run_id;
nFES=popuSize;
popuold=popu;
popold=popu;
fitold=fit;
memory_size = 10;
memory_SF = 0.5 .* ones(memory_size, 1);
memory_CSF = 0.5 .* ones(memory_size, 1);
memory_pos = 1;

while nFES < FES

    mem_rand_index = ceil(memory_size * rand(popuSize, 1));
    mu_SF = memory_SF(mem_rand_index);
    mu_CSF = memory_CSF(mem_rand_index);

    CSF = normrnd(mu_CSF, 0.1);
    term_pos = mu_CSF == -1;
    CSF(term_pos) = 0;
    CSF = min(CSF, 1);
    CSF = max(CSF, 0);
    DSF = ceil(D .* CSF);
    [ti,~] = find(DSF==0);
    DSF(ti) = 1;
    SF = mu_SF + 0.1 * tan(pi * (rand(popuSize, 1) - 0.5));
    pos = find(SF <= 0);
    while ~ isempty(pos)
        SF(pos) = mu_SF(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(SF <= 0);
    end
    SF = min(SF, 1);

    [~,sorted_index] = sort(fit);
    pNP = max(round(0.1 * popuSize), 2);
    randindex = ceil(rand(1, popuSize) .* pNP);
    randindex = max(1, randindex);
    pbest = popu(sorted_index(randindex), :);
    for i = 1:popuSize
        randPopuList = randperm(popuSize);
        randPopuList = setdiff(randPopuList,1,'stable');
        indiR1 = popu(randPopuList(1),:);
        indiR2 = popu(randPopuList(2),:);
        indiR3 = popu(randPopuList(3),:);
        randPopuList = randperm(size(popuold,1));
        randPopuList = setdiff(randPopuList,1,'stable');
        indiR4 = popuold(randPopuList(1),:);
        randDimList = randperm(D);
        randDSFList = randperm(DSF(i));
        seleDim = sort(randDimList(1:randDSFList(1)));
        seleDimLen = length(seleDim);
        if rand > 0.5 * nFES / FES
            if seleDimLen == 1
                s1= SF(i)*HyperSphereTransform_1D(indiR1,indiR2,seleDim);
            elseif seleDimLen == 2
                s1= SF(i)*HyperSphereTransform_2D(indiR1,indiR2,seleDim);
            else
                s1= SF(i)*HyperSphereTransform(indiR1,indiR2,seleDim);
            end
            popu(i,seleDim) = indiR3(seleDim) + s1;
        else
            ss= SF(i)*(pbest(i,:) - popu(i,:) + indiR1 - indiR4);
            s1= ss((seleDim));
            popu(i,seleDim) = popu(i,seleDim) + s1;
        end
    end

    popu = BoundaryDetection(popu, lu);
    fit=Evaluation(popu, benchmark, f);
    fit=fit';
    dif = abs(fitold - fit);
    [fit,II]=min([fitold, fit], [], 2);
    popunew=[popu;popold];
    [~,XI]=unique(popunew,'rows');
    if length(XI)<size(popunew, 1)
        popuold=popunew(XI);
    else
        popuold=popold;
    end
    popu(II==1,:)=popold(II==1,:);
    popold=popu;
    fitold=fit;

    goodCSF = CSF(II == 2);
    goodSF = SF(II == 2);
    dif_val = dif(II == 2);
    num_success_params = numel(goodCSF);
    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        memory_SF(memory_pos) = (dif_val' * (goodSF .^ 2)) / (dif_val' * goodSF);
        if max(goodCSF) == 0 || memory_CSF(memory_pos)  == -1
            memory_CSF(memory_pos)  = -1;
        else
            memory_CSF(memory_pos) = (dif_val' * (goodCSF .^ 2)) / (dif_val' * goodCSF);
        end
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1;
        end
    end

    for i=1:popuSize
        nFES = nFES + 1;
        if fit(i) < bsf_fit_var
            bsf_fit_var = fit(i);
        end
        if mod(nFES, optimalChart_interval) == 0
            optimalChart = [optimalChart;bsf_fit_var];
        else
            if nFES == FES
                optimalChart = [optimalChart;bsf_fit_var];
            end
        end
        if nFES > FES; break; end
    end

    output.popu=popu;
    output.bsf_fit_var=bsf_fit_var;
    output.optimalChart=optimalChart;
    %% Survival
    optimal = bsf_fit_var;
    fprintf('problem %5.0f time %5.0f |%5.0f -----> %9.16f\n',f,run_id,nFES,optimal);
end

function ss = HyperSphereTransform(c,d,pp)
D=length(pp);
A=c(pp)-d(pp);
R=norm(A,2);
O(D-1)=  2*pi*rand ;
for i=1:D-2
    O(i)=  rand*pi  ;
end
C(1)=R*prod(sin(O));
for i=2:D-1
    C(i)= R*cos(O(i-1)) *prod(sin(O(i:D-1)));
end
C(D)=R*cos(O(D-1));
ss=C;
function ss = HyperSphereTransform_1D(c,d,pp)
R=  abs(c(pp)-d(pp));
C = R*cos(2*pi*rand);
ss=C;
function ss = HyperSphereTransform_2D(c,d,pp)
A=c(pp)-d(pp);
R=norm(A,2);
o1=2*pi*rand;
C=zeros(1,2);
C(1) =R*sin(o1);
C(2) =R*cos(o1);
ss=C;