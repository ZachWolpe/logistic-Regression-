
* ________________________________________ LOGISTIC REGRESSION ________________________________________;


data a (keep = age chd )  ;
input nr age chd @@;
cards ;
1 20 0 2  23 0 3  24 0 4  25 0
5 25 1 6  26 0 7  26 0 8  28 0
9 28 0 10 29 0 11 30 0 12 30 0
13 30 0 14 30 0 15 30 0 16 30 1
17 32 0 18 32 0 19 33 0 20 33 0
21 34 0 22 34 0 23 34 1 24 34 0
25 34 0 26 35 0 27 35 0 28 36 0
29 36 1 30 36 0 31 37 0 32 37 1
33 37 0 34 38 0 35 38 0 36 39 0
37 39 1 38 40 0 39 40 1 40 41 0
41 41 0 42 42 0 43 42 0 44 42 0
45 42 1 46 43 0 47 43 0 48 43 1
49 44 0 50 44 0 51 44 1 52 44 1
53 45 0 54 45 1 55 46 0 56 46 1
57 47 0 58 47 0 59 47 1 60 48 0
61 48 1 62 48 1 63 49 0 64 49 0
65 49 1 66 50 0 67 50 1 68 51 0
69 52 0 70 52 1 71 53 1 72 53 1
73 54 1 74 55 0 75 55 1 76 55 1
77 56 1 78 56 1 79 56 1 80 57 0
81 57 0 82 57 1 83 57 1 84 57 1
85 57 1 86 58 0 87 58 1 88 58 1
89 59 1 90 59 1 91 60 0 92 60 1
93 61 1 94 62 1 95 62 1 96 63 1
97 64 0 98 64 1 99 65 1 100 69 1
;
run ;

 
 
proc logistic data=a;
	model chd(event='1') = age / outroc=rocl lackfit;
	output out=n p=pred;
run;



*** ---------------------  Logistic Regression --------------------- ***;
title "Logistic Regression";

proc iml;
use a;
read all into xy;
n=nrow(xy);
x = J(n,1,1)||xy[,1];
y = xy[,2];




* initialize Beta=0 ==> p=1;
bo = {0,0};
diff=999;

start logistic_regression;
	do i=1 to 20 while (diff>0.00001);
		lo = x*bo;
		o = exp(lo);
		p = o/(1+o);
		w = diag(p#(1-p));
		
		bn = bo + inv(x`*w*x)*x`*(y-p);
		print i bo bn;
		diff = abs(bn-bo);
		bo = bn;
	end;
finish logistic_regression;

call logistic_regression;



*** ---------------------  Hosmer-Lemeshow Test --------------------- ***;
title 'Hosmer-Lemeshow Test';

* number of groups;
g=10;

start Hosmer_Lemeshow;
	* get predicted probabilities;
	pred = exp(x*bn) / (1 + exp(x*bn));
	
	d = x || y || pred;
	call sort(d, {4});
	
	
	do i=g to n by g;
		in = i-g+1;
		slice = d[in:i,];
		
		* expected frequency Y=1 per group;
		e = mean(slice[,4])*nrow(slice);
		
		* expected vs actual Y=1;
		ea = ea // (e || sum(slice[,3]));
	end;
	
	* expected vs actual Y=0;
	zeros_ea = (nrow(slice)-ea[,1]) || (nrow(slice)-ea[,2]);
	
	eaa = ea // zeros_ea;
	chisq_stat = 0;
	
	do i=1 to nrow(eaa);
		chisq_stat = 
			chisq_stat + (eaa[i,2]-eaa[i,1])**2 / eaa[i,1];
	end;
	critical_value = quantile('chisq', 0.95, g-2);
	p_value = 1 - cdf('chisq', g-2, chisq_stat);
	
	print chisq_stat critical_value p_value;
finish Hosmer_Lemeshow;

call hosmer_lemeshow;
print ea zeros_ea chisq_stat;



*** ---------------------  ROC Curve --------------------- ***;
title 'ROC Curve';


start confusion_matrix;
	cm = y || classification;
	
	TP=0; TN=0; FP=0; FN=0;
	do c=1 to nrow(cm);
		if cm[c,1] = cm[c,2] then do;
			* TRUE;
			if cm[c,1]=1 then; TP = TP + 1;
			if cm[c,1]=0 then; TN = TN + 1;
		end;
		if cm[c,1] ^= cm[c,2] then do;
			* FALSE;
			if cm[c,2]=1 then; FP = FP + 1;
			if cm[c,2]=0 then; FN = FN + 1;
		end;
	end;
	
	sensitivity = TP/(TP+FN);
	specificity = TN/(TN+FP);
finish confusion_matrix;


start roc;
	pred = exp(x*bn) / (1 + exp(x*bn));
	d = x || y || pred;
	call sort(d, {4});

	do i=1 to 100;
		* for each classification threshold;
		threshold = i/100;
		classification = d[,4] > threshold;
		call confusion_matrix;
		
		roc_data = roc_data // (sensitivity || (1-specificity));
	end;
finish roc;

call roc;

create roc_data from roc_data[colname={'sensitivity' '_1_specificity'}];
append from roc_data;
quit;


proc sgplot data=roc_data;
	series y=sensitivity x=_1_specificity;
	lineparm x=0 y=0 slope=1; 
title "ROC curve";
run;



*** ---------------------  Gini Coefficient --------------------- ***;
title "Gini Coefficient";

*/
ALGORITHM:
	compute line equation for each point {p_i, p_(i+1)} pairing
	randomly sample from domain space
	determine which line segment the point corresponds to
	use Monte Carlo Integration: declare above or below line 
/;


proc iml;
use roc_data;
read all into yx;
y = yx[,1];
x = yx[,2];
n = nrow(yx);
iters = 10000;


start find_equation;
	m = (y2-y1)/(x2-x1);
	c = y2 - m*x2;
finish find_equation;


* determine line equations;
do i=1 to n-1;

	y2 = y[i+1]; y1 = y[i];
	x2 = x[i+1]; x1 = x[i];
	
	* if Y and X differ from Y-1 & X-1;
	if (y2^=y1) & (x2^=x1) then do;
		call find_equation;
		equations  = equations // (x1 || x2 || y1 || y2 || m || c);
	end;
end;




* MC integration;
do i=1 to iters;
	sample_y = rand('uniform');
	sample_x = rand('uniform');
	
	
	* find location w.r.t X;
	loc = 0;
	do i=1 to nrow(equations);
		if sample_x > equations[i,1] then; loc = loc + 1;
	end;
	
	if loc ^= 0 then do;
		* use appropriate equation to evaluate sample;
		fx = equations[loc,5]*sample_x + equations[loc,6];
		if sample_y < fx then; 
			samples = samples // (sample_x || sample_y || 1);
		if sample_y > fx then; 
			samples = samples // (sample_x || sample_y || 0);
	end;
end;

print 'Gini coefficient:' (mean(samples[,3]));



create samples from samples[colname={'x' 'y' 'result'}];
append from samples;
quit;

















