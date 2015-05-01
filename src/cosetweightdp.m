clc; clear; close all;

bin2oct = @(x) str2num(dec2base(bin2dec(x), 8));
oct2int = @(x) base2dec(num2str(x), 8);

CONSTRAINTLENGTH = [3];
CODEGENERATOR = [bin2oct('110') bin2oct('101') bin2oct('111')];
tr_ = poly2trellis(CONSTRAINTLENGTH, CODEGENERATOR);
% feedforward encoders only
tr_.nextStates = tr_.nextStates + 1;
tr_.outputs = arrayfun(oct2int, tr_.outputs);

v = sum(CONSTRAINTLENGTH); k = 1; n = 3; m = max(CONSTRAINTLENGTH) - 1; 
h = 10; T = (m + h); N = n * T; Y = randi([0, 2^n - 1], [1, T]);

tr.activeStates = zeros(T+1, tr_.numStates);
tr.branchMetric = zeros(T, tr_.numStates, tr_.numInputSymbols);

tr.activeStates(1, 1) = 1;
for t = 2:T+1
    if t <= h + 1
        tr.activeStates(t, ...
            unique(tr_.nextStates(tr.activeStates(t-1, :) == 1, :))) = 1;
    else
        tr.activeStates(t, ...
            unique(tr_.nextStates(tr.activeStates(t-1, :) == 1, 1))) = 1;
    end
    tr.branchMetric(t-1, :, :) = ...
        arrayfun(@(x) sum(dec2bin(x) == '1'), bitxor(tr_.outputs, Y(t-1)));
end

dp = zeros(T+1, tr_.numStates, N+1);
dp(T+1, 1, 1) = 1;

for t = T+1:-1:2
    for s = find(tr.activeStates(t-1, :) == 1)
        for u = 1:tr_.numInputSymbols
            if tr.activeStates(t, tr_.nextStates(s, u)) == 1
                for w = 1:N-n+1
                    dp(t-1, s, w+tr.branchMetric(t-1, s, u)) = ...
                    dp(t-1, s, w+tr.branchMetric(t-1, s, u)) + ...
                    dp(t, tr_.nextStates(s, u), w);
                end
            end
        end
    end
end

tr_ = poly2trellis(CONSTRAINTLENGTH, CODEGENERATOR);
CODEWORDS = dec2bin(0:(2^(k*h)-1), k*h) - '0';

cosetwddp = permute(dp(1, 1, :), [1, 3, 2]);
cosetwd = zeros(1, N+1);

Y_ = zeros(1, N);
for i = 1:T
    Y_((i-1)*n+1 : i*n) = dec2bin(Y(i), n) - '0';
end

for i = 1:2^(k*h)
    elem = convenc([CODEWORDS(i, :), zeros(1, m)], tr_);
    wt = sum(bitxor(elem(1:N), Y_));
    cosetwd(wt+1) = cosetwd(wt+1) + 1;
end