# 2007942 - PX457 Assignment 1
# gradient.py

import sys
import math
import scipy.optimize
from typing import Tuple

def f(xy: tuple[float, float]) -> float:
   x, y = xy   # unpack values from tuple
   g = 2*x**3 - 3*x**2 + 5 # expression with x terms from f(x,y) equation
   h = 2*y**3 - 3*y**2 + 5 # expression with y terms from f(x,y) equation
   return g * h
     
def main() -> None:
    # read values from stdin
   x = float(sys.stdin.readline())
   y = float(sys.stdin.readline())
   t = int(sys.stdin.readline())

    # input validation
   if not (0.0 < x < 2.0 and 0.0 < y < 2.0):
        print("Invalid input: x/y outside of range.")
        sys.exit(1)
   if not (1 <= t <= 4):
        print("Invalid input: t outside of range.")
            
    # scipy optimiser
   optima = scipy.optimize.minimize(f, [x,y], method="BFGS")
   if not optima.success:
       print(f"Did not converge within {optima.nit} iterations.")
            
    # unpack results from the minimised functions
   x_minima, y_minima = optima.x
   print(f"Minimum found at {x_minima:.15f} {y_minima:.15f} in {optima.nit} iterations")
         
if __name__ == "__main__":
   main()
