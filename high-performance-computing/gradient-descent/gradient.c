// 2007942 - PX457 Assignment 1 Submission
// gradient.c

#include <stdio.h>
#include <math.h>

int main(void) {
	double x, y; // initial x, y
	int t; // control parameter -- change
	double gamma = 0.01; // learning rate for gradient descent
	double tolerance = 1e-15; // convergence tolerance
	int max_iterations = 200; // maximum iterations for gradient descent

	// read input from stdin, matching
	// the formatting of the assignment
	// specification: x_0 -> y_0 -> t
	scanf("%lf", &x);
	scanf("%lf", &y);
	scanf("%d", &t);
	
	// input validation: 0 < x_0, y_0 < 2.0
	if (x <= 0.0 || x >= 2.0 || y <= 0.0 || y >= 2.0) {
		printf("Invalid input.\n");
		return 1;
	}

	// input validation: 1 < t < 4
	if (t < 1 || t > 4) {
		printf("Invalid t input.\n");
		return 1;
	}
	
	/*
	* 	The gradient descent method uses a recurrence relation
	* x_{n+1} = x_n -  γ*∂f/∂x (1)
	* y_{n+1} = y_n -  γ*∂f/∂y (2)
	* The polynomial takes the form f(x,y) = (2x^3 - 3x^2 + 5)(2y^3 - 3y^2 + 5) (3)
	* We can simplify this to f(x,y) = g(x)h(y)
	* This gives the following partial derivatives:
	* ∂f/∂x = (6x^2 - 6x)(2y^3 - 3y^2 + 5) = g'(x)h(y) (4)
	* ∂f/∂y = (6y^2 - 6y)(2x^3 - 3x^2 + 5) = h'(y)g(x) (5)
	*/
	
	double f, next_f; // values for f(x.y) and f(x_{n+1}, y_{n+1})
	double g, h; // values for g(x), h(y)
	double g_dash, h_dash; // values for g'(x), h'(y)
	double delf_dely, delf_delx; // values for ∂f/∂x, ∂f/∂y
	double next_x, next_y;
	
	// first calculation using formulae above
	g = 2.0*x*x*x - 3.0*x*x + 5.0;
	h = 2.0*y*y*y - 3.0*y*y + 5.0;
	g_dash = 6.0*x*x - 6.0*x;
	h_dash = 6.0*y*y - 6.0*y;
	f = g*h; // Equation (3)
	
	int i; // for loop index
	for (i=1; i <= max_iterations; i++) {
		delf_delx = g_dash*h; // Equation (4)
		delf_dely = h_dash*g; // Equation (5)
		
		// gradient descent using Equations (1) and (2)
		next_x = x - gamma*delf_delx;
		next_y = y - gamma*delf_dely;
		
		// recalculate g(x), h(y) and their derivatives at the new values
		double next_g = 2.0*next_x*next_x*next_x - 3.0*next_x*next_x + 5.0;
		double next_h = 2.0*next_y*next_y*next_y - 3.0*next_y*next_y + 5.0;
		
		double next_g_dash = 6.0*next_x*next_x - 6.0*next_x;
		double next_h_dash = 6.0*next_y*next_y - 6.0*next_y;
		
		next_f = next_g * next_h;
		
		// check for convergence: |delta x|, |delta y|, |delta f| < tolerance
		if (
			fabs(next_x - x) < tolerance &&
			fabs(next_y - y) < tolerance &&
			fabs(next_f - f) < tolerance
		) {
			x = next_x;
			y = next_y;
			f = next_f;
			break;		
		}
		
		// update each variable for the next iteration
		x = next_x;
		y = next_y;
		g = next_g;
		h = next_h;
		g_dash = next_g_dash;
		h_dash = next_h_dash;
		f = next_f;
		
	}
	
	if (i > max_iterations) {
		printf("Error: did not converge within 200 iterations.\n");
		return 1;	
	} else {
		 printf("Minimum found at %.15f %.15f in %d iterations\n", x, y, i);
	}
	
	return 0;
	
}