function index = findridges(Spec,delta)
%Ridge detection algorithm,i.e., Algorithm 1 in paper:Separation of Overlapped Non-Stationary Signals by Ridge Path Regrouping and Intrinsic Chirp Component Decomposition
% ,IEEE Sensors journal,2017.
% The algorithm is originally introduced in the paper: Algorithms for blind components separation and extraction from the time-frequency distribution of their mixture. 
% EURASIP Journal on Advances in Signal Processing, 2004.


%Spec£ºTime-Frequency distribution of the signal
%delta£ºmaximum allowable frequency variation between two consecutive points
%index£ºThe obtained frequency indexs at each time instant

Spec = abs(Spec);
[M,N] = size(Spec);
index = zeros(1,N);
[fmax,tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(tmax) = fmax;

f0 = fmax;

for j = (min(tmax+1,N)):N
    low = max(1,f0-delta);
    up = min(M,f0+delta);
    [~,f0] = max(Spec(low:up,j));
    f0 = f0 + low - 1;
    index(j) = f0;
end

f1 = fmax;
for j = (max(1,tmax-1)):-1:1
    low = max(1,f1-delta);
    up = min(M,f1+delta);
    [~,f1] = max(Spec(low:up,j));
    f1 = f1 + low - 1;
    index(j) = f1;
end
end

