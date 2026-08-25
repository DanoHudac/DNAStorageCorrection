t = [6.5, 181.216, 8924.174];
comb = [582379, 11762970, 537500000];
er = [1, 2, 3];

subplot(2,1,1);
plot(er, t, '-o',"Color",'red');
xlabel('error count');
ylabel('t (m)');
legend('time of running');

subplot(2,1,2);
plot(er, comb,'-o',"Color",'blue');
xlabel('error count');
ylabel('tested combinations');
legend('tested combinations');

sgtitle('k = 32 (1, 2, 3 errors)')
