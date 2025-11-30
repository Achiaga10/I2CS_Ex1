# I2CS_Ex1_proj
I2CS_Ex1_proj is an assignment given by the Ariel University.
The project deals with polynomials.
## Description
This project is used mainly to deal with polynomials in a form of arrays. for example:
{-1.1,2.3,3.1} represents "3.1x^2 +2.3x -1.1" polynomail. 
To deal with each polynomial,This module can handle multiple things like, calculate the f(x),
find the derivative, polynomial multiplication and addition, 
calculate the coefficients of polynomial using 2d points, find area
between two polynomials and more.

## Getting started

## Functions

- _f(double [] polynomial, double x)_ - Computes the f(x) value of the polynomial function at x.
- _root_rec(double[] p, double x1, double x2, double eps)_ - This function computes an x value (x1<=x<=x2) for which |p(x)| < eps.
- _PolynomFromPoints(double[] xx, double[] yy)_ - This function computes a polynomial representation from a set of 2D points on the polynom.
- _equals(double[] p1, double[] p2)_ - returns true if p1 represents the same polynomial function as p2.
- _poly(double[] poly)_ - Computes a String representing the polynomial function.
- _sameValue(double[] p1, double[] p2, double x1, double x2, double eps)_ - This function computes an x value (x1<=x<=x2)
     for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
- _length(double[] p, double x1, double x2, int numberOfSegments)_ - This function computes an approximation of the length of the function between f(x1) and f(x2).
- _area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid)_ - This function computes an approximation of the area between the polynomial functions within the x-range.
- _getPolynomFromString(String p)_ - This function computes the array representation of a polynomial function from a String representation.
- _add(double[] p1, double[] p2)_ - This function computes the polynomial function which is the sum of two polynomial functions (p1,p2).
- _mul(double[] p1, double[] p2)_ - This function computes the polynomial function which is the multiplication of two polynoms (p1,p2).
- _derivative (double[] po)_ - This function computes the derivative of the p0 polynomial function.
- _reverseArray(double[] p)_ - This function reverses a double[] array.


## Usage
```java
package assignments.Ex1;
static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
double[] p = {-1.1,2.3,3.1};
double[] xx = {1, -1, 0};
double[] yy = {3, 3, 2};

// returns 4
double fx1 = Ex1.f(po1, 1);

// returns p1+p2
double[] p12 = Ex1.add(po1, po2);

// returns p1*p2
double[] p12 = Ex1.mul(po1, po2);

// returns derivative po2
double[] po2 = Ex1.derivative(p);

//returns "3.1x^2 +2.3x -1.1"
double[] p1 = Ex1.getPolynomFromString(sp);

//returns true if po2 == p else false
boolean eq = Ex1.equals(po2, p);

//returns intersection point between two polynomials
double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);

// find area between two polynomials in range
double a1 = Ex1.area(po1, po2, x1, x2, 100);

// returns coefficients of the polynomial  
double [] result = Ex1.PolynomFromPoints(xx, yy);

// return length of polynomial in range
double ans = Ex1.length(p, -3.82, 0, 100);
```
## Result
As shown in figure 1, using the above function, It calculated the area between two functions, as well as finding the derivative and more.

![alt text](https://github.com/Achiaga10/I2CS_Ex1/blob/main/src/images/finalresult.png)
Figure 1


## Author
Achia Gabbay
