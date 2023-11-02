clc; clear all;close all;
%%
N = 60; % dim of \rho
phase = 0;%linspace(0,2*pi,365); %squeezing angle
alpha = 5; %coherent amplitude
alpha_minus = -alpha; %coherent amplitude
rdB = 1;
r = rdB/20*log(10); %squeezing factor
%%
a = diag(sqrt(1:N),1); %annihilation  operator
ad = a'; %creation operator

% vacuum state
vacuum = zeros(N+1,N+1);
vacuum(1,1) = 1;

D = expm(alpha*a'-conj(alpha)*a); %displacement operator
D_minus = expm(alpha_minus*a'-conj(alpha_minus)*a); %displacement operator
S = expm(0.5*(r*exp(-1i*phase)*a^2-r*exp(1i*phase)*ad^2)); %squeezing operator

% coherent state
coherent = D*vacuum*D';
coherent_minus = D_minus*vacuum*D_minus';
diag_c = diag(abs(coherent));
Nor_c = trace(coherent);
% figure(1);
% subplot(1,2,1);bar3(abs(coherent))
% subplot(1,2,2);bar(diag_c)
% axis square

% cat state (even)
cat_even = coherent + coherent_minus + D*vacuum*D_minus' + D_minus*vacuum*D';
diag_cat_e = diag(abs(cat_even));
Nor_cat_e = trace(cat_even);
figure(2);
subplot(1,2,1);bar3(abs(cat_even)/Nor_cat_e)
subplot(1,2,2);bar(diag_cat_e/Nor_cat_e)
axis square
title('cat even state (alpha = 1.35)');
% cat state (odd)
cat_odd = coherent + coherent_minus - D*vacuum*D_minus' - D_minus*vacuum*D';
diag_cat_o = diag(abs(cat_odd));
Nor_cat_o = trace(cat_odd);
figure(3);
subplot(1,2,1);bar3(abs(cat_odd)/Nor_cat_o)
subplot(1,2,2);bar(diag_cat_o/Nor_cat_o)
axis square
title('cat odd state (alpha = 1.35)');
%check cat state
Loss = 0;
r = sqrt(1-Loss);
phase_c = pi;
cat = zeros(N,N);
for n = 1 : N+1
    for m = 1 : N+1
        cat(n,m) = (r*alpha)^(n-1)*(r*conj(alpha))^(m-1)/sqrt(factorial(n-1)*factorial(m-1))*exp(-abs(alpha*r)^2)*(1+(-1)^(n+m)+exp(-2*Loss*abs(alpha)^2)*(exp(1i*phase_c)*(-1)^(n-1)+exp(-1i*phase_c)*(-1)^(m-1)))/(2*(1+exp(-2*abs(alpha)^2)*cos(phase_c)));
    end
end
diag_cat = diag(abs(cat));
Nor_cat = trace(cat);
% figure(4);
% subplot(1,2,1);bar3(abs(cat)/Nor_cat)
% subplot(1,2,2);bar(diag_cat/Nor_cat)
% axis square

%%
% single photon
single = ad*vacuum*ad';
% figure;
% subplot(1,2,1);bar3(abs(single))
% subplot(1,2,2);bar(diag(abs(single)))
% axis square

% sq state
% sq = D*S*vacuum*S'*D';
sq = S*vacuum*S';
alpha_s = 5;
D_s = expm(alpha_s*a'-conj(alpha_s)*a); %displacement operator

% sq - 1
sq_a = a*D_s*S*vacuum*S'*D_s'*a';
diag_sq_a = diag(abs(sq_a));
Nor_sq_a = trace(sq_a);
figure(5);
subplot(1,2,1);bar3(abs(sq_a)/Nor_sq_a)
subplot(1,2,2);bar(diag_sq_a/Nor_sq_a)
axis square
title('5dB sq vacuum state - one photon');
% sq + 1
sq_ad = ad*D_s*S*vacuum*S'*D_s'*ad';
diag_sq_ad = diag(abs(sq_ad));
Nor_sq_ad = trace(sq_ad);
figure(6);
subplot(1,2,1);bar3(abs(sq_ad)/Nor_sq_ad)
subplot(1,2,2);bar(diag_sq_ad/Nor_sq_ad)
axis square
title('5dB sq vacuum state + one photon');