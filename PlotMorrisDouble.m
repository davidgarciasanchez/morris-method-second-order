function [] = PlotMorrisDouble(OUT)

maxim=(max(OUT(:,3)))*1.1;
mat=[0:maxim/100:maxim];

p1=md(mat);


x=OUT(:,3) %Mu
y=OUT(:,4) %StdDev

grid on
hold on

axis([0 maxim 0 maxim]);
ylabel('\sigma(h)')
xlabel('\mu (h)')

str1 = num2str(OUT(:,1),'%1d'); % ATTENTION le premier argument de NUM2STR doit être un vecteur
str2 = num2str(OUT(:,2),'%1d'); % ATTENTION le premier argument de NUM2STR doit être un vecteur
str=strcat('w(',str1,',',str2,')');


text(x+0.05,y+0.05,str) 


plot(x,y,'ko',p1,mat,'--b')
print('fig4.eps','-deps', '-FHelvetica:18')

function [mdd] = md(matrix)
    mdd = matrix/sqrt(1);
endfunction
endfunction

