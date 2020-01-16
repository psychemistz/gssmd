% figure 4 & 5
clc
clear all
close all

diff_mean = [0, 1, 3, 5, 10, 20, 30];
nsample = 100;
rng('default')
background = randn(1,nsample);
GSSMD_all = zeros(7, 200);
SSMD_all =  zeros(7, 200);

for i = 1:7
    dm = diff_mean(i);
    
    for j = 1:200
        SNR = 50-0.5*j;
        bgn = awgn(background, SNR);
        tgn = awgn(background + dm,SNR);
        measured_SNR(j)=SNR;
        GSSMD_all(i,j) = gssmd(bgn, tgn);
        SSMD_all(i,j) = ssmd(bgn, tgn);        
    end
end

% plot SSMD
plot(measured_SNR, SSMD_all(1,:))
hold on
plot(measured_SNR, SSMD_all(2,:))
plot(measured_SNR, SSMD_all(3,:))
plot(measured_SNR, SSMD_all(4,:))
plot(measured_SNR, SSMD_all(5,:))
plot(measured_SNR, SSMD_all(6,:))
plot(measured_SNR, SSMD_all(7,:))
hold off
axis([-50 50 -5 30])
%legend("SSMD0", "SSMD1", "SSMD3", "SSMD5", "SSMD10", "SSMD20", "SSMD30", 'Location', 'northwest')

% plot GSSMD
plot(measured_SNR, GSSMD_all(1,:))
hold on
plot(measured_SNR, GSSMD_all(2,:))
plot(measured_SNR, GSSMD_all(3,:))
plot(measured_SNR, GSSMD_all(4,:))
plot(measured_SNR, GSSMD_all(5,:))
plot(measured_SNR, GSSMD_all(6,:))
plot(measured_SNR, GSSMD_all(7,:))
hold off
axis([-50 50 -0.2 1])
%legend("GSSMD0","GSSMD1", "GSSMD3", "GSSMD5", "GSSMD10", "GSSMD20", "GSSMD30", 'Location', 'northwest')

% save results
csvwrite("GSSMD.csv", GSSMD_all);
csvwrite("SSMD.csv", SSMD_all);