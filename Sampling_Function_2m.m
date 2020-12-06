function [Outmatrix, OutFact, OutPlace] = Sampling_Function_2m(p, k, r, UB, LB, GroupMat)

#Modification by David Garcia Sanchez 2010. 

# The New Morris Method was proposed by Campolongo and Braddock [Reliab. Engng Syst. Saf. 64 (1999) 1] as an extension of the Morris Method [Technometrics 33 (1991) 161] 
# to include estimation of two-factor interaction effects. 

# Parameters and initialisation of the output matrix

runs=r;
sizea = k;
Delta = p/(2*p-2);
#Delta = 1/5;
NumFact = sizea;
GroupNumber = size(GroupMat,2);

if GroupNumber ~ 0;
    sizea = size(GroupMat,2);
endif

sizeb = sizea + 1;
sizec = sizea - 1;
Outmatrix = [];
OutFact = [];
OutPlace =[];

# For each i generate a trajectory  
for i=1:r
    for l=1:sizea/2;
	
		randmult = ones(k,1);           
    
		v = rand(k,1);                  
		randmult (find(v < 0.5))=-1;
		randmult = repmat(randmult,1,k);
		DD0 = randmult .* eye(k);     

		# construct the multiple trayectoires for b
    
		H = handcuffed(sizea,l)
		B = createM(H)
    
		# Construct A0, A
		A0 = ones(sizeb,1);
		A = ones(sizeb,NumFact);   

		P0= diag(ones(1,sizea));   

		if GroupNumber ~ 0
			B = B * (GroupMat*P0')';
		endif   

		AuxMat = Delta*0.5*((2*B - A) * DD0 + A);
    
		# a --> Define the random vector x0 for the factors. Note that x0 takes value in the hypercube
		# [0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta] 
		MyInt = repmat([0:(1/(p-1)):(1-Delta)],NumFact,1);     # Construct all possible values of the factors               
   
		v = repmat(rand(NumFact,1),1,size(MyInt,2)+1);     
		IntUsed = repmat([0:1/size(MyInt,2):1],NumFact,1); 
		DiffAuxVec = IntUsed - v;                          
    
		for ii = 1:size(DiffAuxVec,1)
			w(1,ii) = max(find(DiffAuxVec(ii,:)<0));       
		endfor
		x0 = MyInt(1,w)';                                 

		if GroupNumber ~ 0
			B0 = (A0*x0' + AuxMat);
		else
			B0 = (A0*x0' + AuxMat)*P0;
		endif
    
		# c --> Compute values in the original intervals
		# B0 has values x(i,j) in [0, 1/(p -1), 2/(p -1), ... , 1].
		# To obtain values in the original intervals [LB, UB] we compute LB(j) + x(i,j)*(UB(j)-LB(j))
		In = repmat(LB,1,sizeb)' + B0 .* repmat((UB-LB),1,sizeb)';
    
		#  Create the (k-1)-by-k matrix of EEj necessary to compute using the actual Handcuffed vector H(1,sizea) and the results matrix
    
		EEj=zeros(sizea-1,sizea);
    
		for i=1:sizea-1
			#create a matrix EEj of (k-1)-by-k using only the first k-1 results of In (wich is the B* matrix with the delta factors)
			EEj(i,:)= In(i,:);
	
			# take the exact value should change from In(B*) and place in the j change in EEj
			EEj(i,H(i+1))=In(2+i,H(i+1));
   
			#we construct also vectors who change for use in Fact variable down
			Factj(i)=H(i+1);
    
		endfor
    
		Fact = [0 H(1:sizea)];
   
		OutPlace1 = place(sizea,1);
		OutPlace = [OutPlace; OutPlace1];
		Outmatrix = [Outmatrix; In; EEj];
		OutFact = [OutFact; Fact'; Factj'];
    endfor
endfor

endfunction

# functions added by David Garcia 2010 

function [B] = createM(A)

# function that create morris block matrix given a trayectory this should be able to compute EEij and EEi Number of factors
k=size(A,2);
tam=size(A,2);
tamb=tam+1;

B=zeros(k+1,k);

B(1:k+1>=3,A(1))=1;
B(1:k+1>=3,A(2))=1;
B(1:k+1>=4,A(3))=1;
B(1:k+1>=5,A(4))=1;

# El paso del conejo
B(2,A(1))=1;

	for i=6:tamb
		B(1:tamb>=i,A(i-1))=1;
	endfor
	
endfunction

function [a] = handcuffed(sizea, n)

# Hand-cuffed series a series of blocks are obtained such as a <...s+2,2s-1,s+1,0,s,1,s-1,2,s-2...> are obtained
# Mendelshon NS. Handcuffed designs. discrete mathematics 1977.  NOTE: that sizea = k (number of factor) and not s=k/2

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
function [P_matrix] = place(k,R)

P_matrix = [];
for i=1:R
    for j=1:k+1
        P1(j)=j;
    endfor
    
    for j=1:k-1
        P2(j)=88;
    endfor
    
	P_matrix = [P_matrix; P1'; P2'];
endfor

endfunction

