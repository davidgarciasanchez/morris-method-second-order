# morris-method-second-order

### Description
The New Morris Method was proposed by Campolongo and Braddock [Reliab. Engng Syst. Saf. 64 (1999) 1] as an extension of the Morris Method [Technometrics 33 (1991) 161] to include estimation of two-factor interaction effects.

I made this code to sample and calculate the sensitivitiy analysis with two factors interaction effects using a solution of handcuffed prisoners proposed by Campolongo and based on the work of Mendelshon NS. Handcuffed designs. discrete mathematics 1977.

Code and/or part of this code has been used in the following scientific publications:

- Garcia Sanchez D, Lacarrière B, Musy M, Bourges B. Application of sensitivity analysis in building energy simulations: combining first-and second-order elementary effects methods. Energy and Buildings, 2014;68:741–50. (https://www.sciencedirect.com/science/article/abs/pii/S0378778812004768)
- Awad M, Senga Kiesse T, Assaghir Z, Ventura A. Convergence of sensitivity analysis methods for evaluating combined influences of model inputs, Reliability Engineering & System Safety, Volume 189, September 2019, Pages 109-122 (https://www.sciencedirect.com/science/article/abs/pii/S0951832018305763)

### Use

1. [R1 R2] = Sampling_Function_2m("$Pnumber","$Knumber", "$Rnumber", UM(:,1),UM(:,2),ones("$Knumber",0));
2. Run your code/model/experiment based on the trajectoires on R1 output and recover the matrix of results (Res)
3. [OUT] = Morris_measure_double("$Knumber", R1, Res,"$Pnumber"); 
4. Use PlotMorrisDouble to plot or use your plot code on the results matrix [OUR]
5. plotMorrisDouble2(OUT);
### References
https://www.sciencedirect.com/science/article/abs/pii/S0951832002001096

