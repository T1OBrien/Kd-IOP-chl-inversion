function [bbp, bb] = MM01bbpchl(lambda, chl)
% MM01 bbp from chlorophyll
% chl in ug/l or mg/m3
% lambda in nm
% Terence O'Brien
% version 2026/10/16
% follows Morel & Maritorena 2001

N = length(lambda);

for i = 1:length(chl)

bp_550_chl(i) = 0.416.*chl(i).^0.766;

if chl(i) >= 0.02 && chl(i) <= 2 %
    nu(i) = 1/2 * (log10(chl(i))-0.3);
    bbp(i,1:N) = (0.002 + 0.01*(0.50 - 0.25.*log10(chl(i))).*(lambda./550).^nu(i)) .* bp_550_chl(i);
else
    nu(i) = 0;
    bbp(i,1:N) = (0.002 + 0.01*(0.50 - 0.25.*log10(chl(i))).*(lambda./550).^nu(i)) .* bp_550_chl(i);
else 
    disp('chl too low!')
    bbp = 0;
end



Temp = 10; 
Salt = 33; % estimates

[~,bbw]=IOP_of_sw(lambda,Temp,Salt); % IOPs of sw from EB

bb = (1/2).*bbw + bbp;
end
