# I2CS_Ex1_proj
I2CS_Ex1_proj is an assignment given by the Ariel University.
The project deals with polynomials.
## Description
This project is used mainly to deal with polynomials, calculate the f(x),
find the derivative, polynomial multiplication and addition, 
calculate the coefficients of polynomial using 2d points, find area
between two polynomials and more.

## Getting started

### Dependencies

- java openjdk 21.0.5 2024-10-15 LTS
  OpenJDK Runtime Environment Temurin-21.0.5+11 (build 21.0.5+11-LTS)
  OpenJDK 64-Bit Server VM Temurin-21.0.5+11 (build 21.0.5+11-LTS, mixed mode, sharing)

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

## Author
Achia Gabbay