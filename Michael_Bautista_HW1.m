%%HW#1 Name: Michael Bautista, email:mbautis000@citymail.cuny.edu
%%Question 1
clf
m = 317;
a = 151;
r = ones(1,m-1); %create array of temporary values from 1 to 316 

for i=2:m-1
    r(i) = mod(a*r(i-1), m); %replace temp. values using the formula
end

N = 1:m-1;
k = zeros(1,m-1); %create 1 by 316 matrix of 0's

for i=1:m-1
   k(i) = nnz(r==i); %nnz returns non-zero elements for r equals i
end

plot(N,k,'-', 'color', 'red')
title('Relationship between K and N(k)')
xlabel('k')
ylabel('N(k)')
%All N(k) = 1 because each value of k only shows up once 

%%Question 2
% 1/n * summation_from_1_to_n(Uk) for 100 values
% PRNG used in q1:
close all; clear all; clc
m = 317;
a = 151;
% mod(163*192, 347) = 66, which shows that it's full period
r = ones(1,m-1);
for i=2:m-1
    r(i) = mod(a*r(i-1), m);
end
q = r/m;
sum(q)/m
%other PRNG:
method1 = RandStream('mcg16807');
U1 = rand(method1,1,10000)
method2 = RandStream('swb2712');
U2 = rand(method2,1,10000)
method3 = RandStream('mt19937ar');
U3 = rand(method3,1,10000)

% %subplot(1,3,1); semilogx(q)
% subplot(1,3,1); semilogx(U1)
% subplot(1,3,2); semilogx(U2)
% subplot(1,3,3); semilogx(U3)

subplot(2,2,1); semilogx(q); hold on
title('Method from Q1')

subplot(2,2,2);
semilogx(U1); hold on; title('Method 1')
subplot(2,2,3); semilogx(U2); title('Method 2')
subplot(2,2,4); semilogx(U3); title('Method 3')


%The graphs converge to E[U^n] = 1/(n+1), so E[U] = 1/2, and E[U^5] = 1/6

%% Question 3
%Using code from Q2:
close all; clear all; clc
m = 317;
a = 151;
r = ones(1,m-1);
for i=2:m-1
    r(i) = mod(a*r(i-1), m);
end
q = r/m;
sum(q)/m;
method1 = RandStream('mcg16807');
U1 = rand(method1,1,10000);
method2 = RandStream('swb2712');
U2 = rand(method2,1,10000);
method3 = RandStream('mt19937ar');
U3 = rand(method3,1,10000);

correlationcoefficient(U1,U1)
correlationcoefficient(U2,U2)
correlationcoefficient(U3,U3)

subplot(1,3,1); plot(U1(1:2:1999),U1(2:2:2000),'.');
title('Method1')
subplot(1,3,2); plot(U2(1:2:1999),U2(2:2:2000),'.');
title('Method2')
subplot(1,3,3); plot(U3(1:2:1999),U3(2:2:2000),'.');
title('Method3')

function [r] = correlationcoefficient(firstvalue,secondvalue)
    x = firstvalue(1:2:1999);
    y=  secondvalue(2:2:2000);
    b = mean(x);
    c = mean(y);
    b1 = x-b;
    c1 = y-c;
    bcproduct = b1.*c1;
    b1square = b1.^2;
    c1square = c1.^2;
    sumproduct = sum(bcproduct);
    sumbsquare = sum(b1square);
    sumcsquare = sum(c1square);
    product = sumbsquare * sumcsquare;
    denominator = sqrt(product);
    r = sumproduct / denominator
end

%%Question 4
%Using code from Q2:
close all; clear all; clc
m = 317;
a = 151;
r = ones(1,m-1);
for i=2:m-1
    r(i) = mod(a*r(i-1), m);
end
q = r/m;
sum(q)/m;
% other PRNG:
method1 = RandStream('mcg16807');
U1 = rand(method1,1,10000);
method2 = RandStream('swb2712');
U2 = rand(method2,1,10000);
method3 = RandStream('mt19937ar');
U3 = rand(method3,1,10000);

%Looking through columns to satisfy: V1 < 1/4, 1/4 < V16 < 1/2, and 1/2 < V28
V1 = []; V16 = []; V28 = []; 

CubeA = rand(method1,30,30);

for k=1:30
    %Look through row 1 to see if values are < .25
    if CubeA(1,k) < .25
        V1(k) = true;
    else
        V1(k) = false;
    end
    
    %Look through row 16 to see if values btwn .25 and .5
    if CubeA(16,k) > .25 & CubeA(16,k) < .5
        V16(k) = true;
    else
        V16(k) = false;
    end
    
     %Look through row 28 to see if values are > .5
    if CubeA(28,k) > .5
        V28(k) = true;
    else
        V28(k) = false;
    end    
end

total = sum(V1) + sum(V16) + sum(V28);
PrbE = (total./90) %Prb of Method A = .311

%Repeating the process for CubeA = rand(method2,30,30) results in PrbE = .311
%and CubeA = rand(method3,30,30); = .322. The product of all 3 Prb is .0312

%%4b
Cube30for1 = rand(method1,30,10000); %30 dim. cube that uses method 1
Cube30for2 = rand(method2,30,10000);
Cube30for3 = rand(method3,30,10000);

%Checking the # of pseudo random vectors in set E
count1 = 0; count2 = 0; count3 = 0;
for col = 1:10000
    if Cube30for1(1,col) <.25 & Cube30for1(16,col) > .25 & Cube30for1(16,col) < .5 & Cube30for1(28,col) > .5
        count1 = count1 + 1;
    end
    if Cube30for2(1,col) <.25 & Cube30for2(16,col) > .25 & Cube30for2(16,col) < .5 & Cube30for2(28,col) > .5
        count2 = count2 + 1;
    end
    if Cube30for3(1,col) <.25 & Cube30for3(16,col) > .25 & Cube30for3(16,col) < .5 & Cube30for3(28,col) > .5
        count3 = count3 + 1;
    end
end

fraction1 = count1 / 10000 %Fraction of which vectors are in V for method 1
fraction2 = count2 / 10000
fraction3 = count3 / 10000

%fraction1 = .0324, fraction2 = 0, fraction3 = .0332. These 3 probabilities converge to .0321 from (a)

%%4c
%[h,p] = kstest(U1)
 test_cdf = makedist('tlocationscale','mu',75,'sigma',10,'nu',1);
 [h,p] = kstest(fraction3,'CDF',test_cdf,'Alpha',0.1)
 %The returned value of h is 1, which tells us that the kstest rejects the null hypothesis
 %at the 0.0844 level (p = 0.0844) for method 1 (p is the same value for method 2 and 3 as well)
