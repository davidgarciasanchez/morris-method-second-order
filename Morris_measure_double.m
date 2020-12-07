#octave

function [OutMatrix] = Morris_measure_double(NumFact, Sample, Output, p);

Delt = p/(2*p-2)
sizea=NumFact;

#The number of single trayectoires in MT required to compute a complete EEij
#Analysis is given by the formula r=k/2

totalT=size(Sample,1)/((2*sizea)*(sizea/2));

r=sizea/2;

# We use g as a index total index of runs.
g=0;

#we start a loop for the TOTAL number of trayectoires (shared as RBIG)
for i=1:totalT
    
    # We take a single MT (k+1)+(k-1)-by-k matrix of sample values and
    # their respectives results

    stepa=((2*sizea)*(sizea/2))-2;
    stepb=((2*sizea)*(sizea/2))-1;
    
    aa = i+(i-1)*(stepa+1);
    bb = i+(i-1)*stepb+stepb;
    
    Single_Sample = Sample(aa:bb,:);
    Single_OutValues  = Output(aa:bb,:);    
        
    for j=1:sizea/2; #loop for k/2 number of B* necessary to compute EEij
        
        sizemt=(2*sizea)-1;
  
        A = Single_Sample(j+(j-1)*sizemt:j+(j-1)*sizemt+sizemt,:);
        R = Single_OutValues(j+(j-1)*sizemt:j+(j-1)*sizemt+sizemt,:);
        

        H = handcuffed(sizea,j);
        
        Change_factor = zeros(sizea-1,2);
        
        for k=1:sizea-1
            
            if H(k)<H(k+1)
                Change_factor(k,1:2) = [H(k) H(k+1)] ;
            else
                Change_factor(k,1:2) = [H(k+1) H(k)] ;
            endif
            
            g=g+1;
      
            SAmeas(g,1)=Change_factor(k,1);
            SAmeas(g,2)= Change_factor(k,2);
            
            #debug
            SAmeasd(g,1)=Change_factor(k,1);
            SAmeasd(g,2)= Change_factor(k,2);
            
             EE1=A(k,:);
             EE2=A(k+1,:); #i matrix
             EE3=A(k+(sizea+1),:); #j matrix
           
             
             VER1(g,:)=EE1-EE2;
             VER2(g,:)=EE1-EE3;  
            
            
            EEi=(R(k+1)-R(k))/Delt;
            EEj=(R(k+(sizea+1))-R(k))/Delt;
            
            
            SEEij=(R(k+2)-R(k));
            SAmeasd(g,3:5)= [SEEij EEi EEj];
            
            SAmeas(g,3) = abs((R(k+2) - R(k+1) - R(k+(sizea+1)) +R(k)) /(Delt*Delt));
        endfor
        
    endfor
    
endfor

VER1;
VER2;
SAmeas;
SAmeasd;


Nint = (sizea*(sizea-1))/2;
SAmeasN = sortrows(SAmeas,[1 2]);

Tbass=totalT-1;

for y=1:Nint
a = y+(y-1)*Tbass;
b = y+((y-1)*Tbass)+Tbass;
Mut(y) = sum(abs(SAmeasN(a:b,3)))/totalT;
StDevt(y) = ((sum((SAmeasN(a:b,3) - Mut(y)).^2))*(1/(totalT-1))).^0.5;
endfor

Mu=Mut'
StDev=StDevt';
Factors1=SAmeas(1:Nint,1:2);
Factors=sortrows(Factors1);

Factors

OutMatrix=([Factors, Mu, StDev]);

endfunction


function [a] = handcuffed(sizea, n)
# Hand-cuffed series a series of blocks are obtained such as
# a <...s+2,2s-1,s+1,0,s,1,s-1,2,s-2...> are obtained
# Mendelshon NS. Handcuffed designs. discrete mathematics 1977.
# NOTE: that sizea = k (number of factor) and not s=k/2

sf=sizea;
j=1;
k=0;
kk=0;
for s=1:sf;
    if mod(s,2) == 0
        #number if even
        a(s)=a(s-1)-1;
    else
        #number is odd
        a(s)=(sizea-1)-j;
        j=j+1;
    endif
endfor

for s=1:sf
    if mod(s,2)==0
        #even
        a(sizea-k)=(a(sf))-kk;
        k=k+2;
        kk=kk+1;
    endif
endfor
a(2)=sizea-1;

# This indicate the rotation of handcuffed if is necessary

if n>1
    for t=1:n-1
        for j=1:sizea
            a(j)=a(j)+1;
        endfor
         [C,I]=max(a);
         a(I)=C-sizea;
    endfor
endif

#this add the real values to matrix a

for i=1:k
    a(i)=a(i)+1;
endfor

endfunction
    
