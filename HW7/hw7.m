% 1. Calculate the excess returns of conducting carry trade for each of the
% currency pairs

CHF = CHF1MTD156N(2610:8155,2);
EUR = EUR1MTD156N(:,2);
GBP = GBP1MTD156N(3393:8938,2);
JPY = JPY1MTD156N(3393:8938,2);
USD = USD1MTD156N(3393:8938,2);
MTDmat = [CHF EUR GBP JPY USD];
naIndex = ismissing(MTDmat);

% replace missing values in data
for i = 1:length(naIndex)
    if (naIndex(i) == 1)
        MTDmat(i,:) = MTDmat(i-1,:);
    end
end

% calculate carry trade for each of the currency pair
dimensions = size(MTDmat);
row_num = dimensions(1);
col_num = dimensions(2);
ers = zeros(row_num, col_num);


for r =  1:row_num
    index = 1;
    for i = 1:col_num
        for j = i+1:col_num
            if (i ~= j) && (MTDmat{r,i} > MTDmat{r,j})
                ers(r,index) = MTDmat{r,i}  - MTDmat{r,j};
                index = index + 1;
            end
            if (i ~= j) && (MTDmat{r,i} < MTDmat{r,j})
                ers(r,index) = MTDmat{r,j}  - MTDmat{r,i};
                index = index + 1;
            end
        end
    end
end


% Test if the excess returns are significant

s_mean = mean(ers);
s_std = std(ers);
t  = s_mean./(s_std/sqrt(5546));

% t =
% 
%   107.4444  113.1476   76.9714   84.5528  102.3460   81.6712  102.4132   85.4463   81.3323   74.8183
% The t-values showed that the  excess returns are significantly different from zeros.

% 2. Test if beta is 1 for 9 currency pairs

% US JPY
usJPY = 1./DEXJPUS{3914:12851,2};
usJPY = removeNA(usJPY);
rUSD = USD1MTD156N{:,2};
rUSD = removeNA(rUSD);
rJPY = JPY1MTD156N{:,2};
rJPY = removeNA(rJPY);
sJPY = delta_s(usJPY);

x = (rUSD - rJPY)/1200;
x = x(1:end-22);
x = [ones(length(x),1) x];

paramsJPY = inv(x.'*x)*x.'*sJPY;
% [0.003496806847093;-1.106852267123259]

usSZ = 1./DEXSZUS{1:1302,2};
usSZ = removeNA(usSZ);
rUSD = USD1MTD156N{7637:end,2};
rUSD = removeNA(rUSD);
rCHF = CHF1MTD156N{6854:8155,2};
rCHF = removeNA(rCHF);
sSZ = delta_s(usSZ);

x = (rUSD - rCHF)/1200;
x = x(1:end-22);
x = [ones(length(x),1) x];

paramsSZ = inv(x.'*x)*x.'*sSZ;
%[-0.002949627110183;1.555340667365140]
% US EU
usEU = DEXUSEU{1:5546,2};
rEUR = EUR{:,1};
rUSD = USD{:,1};
usEU = removeNA(usEU);
rEUR = removeNA(rEUR);
rUSD = removeNA(rUSD);
sEU = delta_s(usEU);
x = (rUSD - rEUR)/1200;
x = x(1:end-22);
x = [ones(length(x),1) x];

paramsEU = inv(x.'*x)*x.'*sEU;
% [3.831241854915842e-04;-1.441642202670693]
% US UK
usUK = DEXUSUK{3914:12851,2};
rUK = GBP1MTD156N{:,2};
rUSD = USD1MTD156N{:,2};
usUK = removeNA(usUK);
rUK = removeNA(rUK);
rUSD = removeNA(rUSD);
sUK = delta_s(usUK);
x = (rUSD - rUK)/1200;
x = x(1:end-22);
x = [ones(length(x),1) x];

paramsUK = inv(x.'*x)*x.'*sUK;
% [-7.125562887990108e-04;-0.274948797783444]

% the calculated parameters showed that the theory does not hold.
% 3.  Apply trading strategy and test if excess returns are significant.
% in-sample data
erUK = zeros(8635,1);
for i = 1:8635
    x_s = x(1:259+i,:);
    y = sUK(1:259+i);
    size(x_s)
    size(y)
    params = inv(x_s.'*x_s)*x_s.'*y;
    s_exp = [1 x(261)]*params;
    if (s_exp > x(259+i))
        erUK(i) = s_exp - x(259+i);
    else 
        erUK(i) = x(259+i) - s_exp;
    end
end
    
s_mean = mean(erUK);
s_std = std(erUK);
t_UK  = s_mean./(s_std/sqrt(8635));

% t-UK is 72.96, which mean the excess profit is  significant. 

% define a function, replace missing values
function naremoved = removeNA(x)
    naIndex = ismissing(x);

    % replace missing values in data
    for i = 1:length(naIndex)
        if (naIndex(i) == 1)
            x(i,:) = x(i-1,:);
        end
    end
    naremoved = x
end
    

% define a function that calculates s{t+k} -  s{t}
function delta_s = delta_s(x)
    delta_s = zeros(length(x)-22,1);
    for i = 1:length(x)-22
        delta_s(i) = log(x(i+22)) - log(x(i));
    end
end