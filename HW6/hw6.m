% Compare long, short and arbitrage portfolios
% rank portfolios based on culmulative return during past 6 and 9 months
% long 3 top portfolios and short 3 bottom portfolios

J = [6; 9]

IndustryPortfolios = IndustryPortfolios{:,:} * 0.01 % return in %
% calcute culmulative returns

acc_rets = zeros(1117,30)
for i = 8:1124
    rets = ones(1,30);
    for j = 2:7
        rets = rets.*(1+IndustryPortfolios{i-j,:});
    end
    acc_rets(i-7,:) = rets;
end


% go through the sample, add three columns to original dataset: long_ret,
% short_ret
    
        
        