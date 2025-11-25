package assignments.Ex1;

import java.util.Random;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};

	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0; i < poly.length; i++) {
			double c = Math.pow(x, i);// x in the power of i
			ans += c*poly[i]; // calculates the final ans
		}
		return ans;
	}

    /** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0. Basically finds when the function goes from (+) value to (-) and vise versa
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}

    /**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
     * The function calculates the coefficients using the Gauss algorithm for matrix.
	 * @param xx - array of x's values that goes accordingly with yy, for exp: {1,-6,0.3}
	 * @param yy - array of y's values that goes accordingly with xx, for exp: {0,2,-9}
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null; // {x,x2,x3}
		int lx = xx.length;// {y,y2,y3}
		int ly = yy.length;
        double [][] matrix = new double[lx][lx+1];
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            for (int i = 0; i < lx; i++) {
                if (lx == 3){
                    matrix[i] = new double[]{Math.pow(xx[i], 2), xx[i], 1, yy[i]};
                } else{
                    matrix[i] = new double[]{xx[i], 1, yy[i]};
                }
            }
            for (int i=0; i < matrix.length - 1; i++) {
                if(matrix[i][i]==0){
                    for (int j=i+1; j < matrix.length; j++) {
                        if(matrix[j][i]!=0) {
                            double[] tmp = matrix[i];
                            matrix[i] = matrix[j];
                            matrix[j] = tmp;
                            break;
                        }
                    }
                }

                double divisor = matrix[i][i];
                for (int z=0; z < matrix[i].length; z++) {
                    matrix[i][z] = matrix[i][z] / divisor;
                }
                for (int j = 1+i; j < matrix.length; j++) {
                    double makeItZero = matrix[j][i] / matrix[i][i];
                    for(int k=i; k < matrix[j].length; k++){
                        matrix[j][k] = matrix[j][k] - (makeItZero * matrix[i][k]);
                    }
                }
            }
            // Since there are only up to three coefficients,
            // It calculates it manually and insert it into ans.
            ans = new double[lx];
            if (lx == 3){
                ans[2] = matrix[2][3] / matrix[2][2];
                ans[1] = matrix[1][3] - ans[2]*matrix[1][2];
                ans[0] = matrix[0][3] -(ans[1]*matrix[0][1] + ans[2]*matrix[0][2]);
            }
            else {
                ans[1] = matrix[1][2] / matrix[1][1];
                ans[0] = matrix[0][2] - ans[1]*matrix[0][1];
            }

		}
        //It reverses the array to represent in the correct form
		return reverseArray(ans);
	}

    /** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true if p1 represents the same polynomial function as p2.
     */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        if(p1==null && p2==null) return true;

        Random rand = new Random();
        int n = Math.max(p1.length, p2.length) + 1;
        int[] testArray = new int[n];

        for(int i=0; i<n; i++){
            testArray[i] = rand.nextInt(201) - 100;
        }
        for (int i = 0; i < testArray.length; i++) {
            double resultP1 = f(p1,testArray[i]);
            double resultP2 = f(p2,testArray[i]);

            if (Math.abs(resultP1 - resultP2) > EPS) {
                ans = false;
            }

        }
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            ans = ans + poly[0];
            if (poly[0] > 0) {ans ="+" + poly[0];}
            for (int i = 1; i < poly.length; i++) {
                if (i == 1){
                    ans = poly[i]+ "x " + ans;
                    if (poly[i] >= 0 && i < poly.length-1) {
                        ans = "+" + ans;
                    }
                }
                else{
                    ans =  poly[i] + "x^"+ (i) + " " + ans;
                    if (poly[i] >= 0 && i < poly.length-1) {
                        ans = "+" + ans;
                    }
                }
            }
		}
		return ans;
	}

    /**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double increment = 0.000001;
        for (double i = x1; i < x2; i+=increment) {
            if (Math.abs(f(p1,i) - f(p2,i)) < eps){
                return i;
            }
        }
		return x1;
	}

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        if (numberOfSegments <= 0 || x1 >= x2) {
            if (x1 == x2) return 0.0;
            return 0.0;
        }

        double deltaX = (x2 - x1) / numberOfSegments;
        double ans = 0.0;
        double currentY = f(p, x1);

        for (int i = 1; i <= numberOfSegments; i++) {
            double nextX = x1 + i * deltaX;
            double nextY = f(p, nextX);
            double deltaY = nextY - currentY;

            //Math.sqrt(deltaX * deltaX + deltaY * deltaY)
            double segmentLength = Math.hypot(deltaX, deltaY);
            ans += segmentLength;

            currentY = nextY;
        }
        return ans;
    }

	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * And for extra accurecey when there is an intersection point it computes two triangles from both side
     * of the intersection.
     * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
    public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        if (numberOfTrapezoid <= 0 || x1 == x2) return 0.0;
        double h = Math.abs(x2 - x1) /  numberOfTrapezoid;
        if (Math.abs(f(p1,x1) - f(p2,x1)) < EPS) {
            double startBase = Math.abs(f(p1,x1+h) - f(p2,x1+h));
            ans = startBase*h/2;
            System.out.println("Starting with crossed point");
            x1 += h;
        }

        double middleX;
        for (double i = x1; i < x2; i+=h) {
            middleX = sameValue(p1, p2, i, i+h, 0.00001);
            if (middleX != i){
                double baseA = Math.abs(f(p1,i) - f(p2,i));
                double ha = Math.abs(middleX - i);
                ans += baseA * ha / 2;

                double baseB = Math.abs(f(p1,i+h) - f(p2,i+h));
                double hb = Math.abs(i+h-middleX);
                ans += baseB * hb / 2;
            }
            else{
                double a =  Math.abs(f(p1,i) - f(p2,i));
                double b =  Math.abs(f(p1,i+h) - f(p2,i+h));

                ans += (a+b) * h / 2;
            }

        }
        return ans;
    }

    /**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return double [] representing the polynomial function: for: '-1.0x^2 +3.0x +2.0' -> {2.0, 3.0, -1.0}
	 */
    public static double[] getPolynomFromString(String p) {
        String[] words = p.split("\\s+");
        double [] ans = new double[words.length];//  -1.0x^2 +3.0x +2.0 {2.0, 3.0, -1.0}
        String charactersToRemove = "[x^+]";
        for (int i = words.length - 1; i >= 0; i--) {
            words[i] = words[i].split("x")[0];
        }
        for(int i = 0; i < words.length / 2; i++)
        {
            String temp = words[i];
            words[i] = words[words.length - i - 1];
            words[words.length - i - 1] = temp;
        }
        for (int i = 0; i < words.length; i++) {
            ans[i] = Double.parseDouble(words[i]);
        }

        return ans;
    }

	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     * It does it by iterating the arrays and adding the matching index values into one array.
	 * @param p1 first polynomial
	 * @param p2 second polynomial
	 * @return sum of two polynomial functions (p1,p2)
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        int len = Math.max(p1.length, p2.length);
        if (len == 0) return ans;
        ans = new double[len];
        for (int i = 0; i < len; i++) {
            if (i > p2.length - 1) {
                ans[i] = p1[i];
            } else if (i > p1.length - 1) {
                ans[i] = p2[i];
            }
            else{
                ans[i] = p1[i]+p2[i];
            }
        }
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     * It creates a matrix that is p1.length*p2.length+p1.length-1, then multiply each value at p1 with
     * each value of p2 then stores the output in the matrix. Then it adds all matrix values to get the
     * final array of polynomial
     * @param p1 first polynomial
     * @param p2 second polynomial
	 * @return multiplication of two polynoms (p1,p2)
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
		double [][] polysToAdd = new double[p1.length][p2.length+p1.length-1];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                polysToAdd[i][i+j] = p1[i]*p2[j];
            }
        }
        ans = add(ans, polysToAdd[0]);
        for (int i = 1; i < p1.length; i++) {
            ans = add(ans, polysToAdd[i]);
        }
		return ans;
	}

    /**
	 * This function computes the derivative of the p0 polynomial function.
     * It iterates po from the po[1] index, then each value is multiplied by i which represents the power.
	 * @param po
	 * @return array of the derivative of po
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;
        if (po.length == 0) return ans;
        ans = new double[po.length - 1];
        for (int i = 0; i < ans.length; i++) {
            ans[i] = po[i+1]*(double)(i+1);
        }
		return ans;
	}

    /**
     * This function reverses a double[] array
     * @param p double[]
     * @return p reversed double[]p
     */
    public static double[] reverseArray(double[] p) {
        for (int i = 0; i < p.length / 2; i++) {
            double temp = p[i];
            p[i] = p[p.length - i - 1];
            p[p.length - i - 1] = temp;
        }
        return p;
    }
}
//example of poly as array and string
// double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
