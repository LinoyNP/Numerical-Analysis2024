def my_secant_method(f, TOL=1e-6, N=50):
    x0=20
    x1=40
    #print("{:<10} {:<15} {:<15} {:<15}".format("Iteration", "xo", "x1", "p"))
    for i in range(1,N+1):
        if f(x1) - f(x0) == 0:
            print( " method cannot continue.")
            return  None

        p = x0 - ((x1 - x0) / (f(x1) - f(x0))) * f(x0)

        if abs(p - x1) < TOL:
            return p  # Procedure completed successfully
        #print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f}".format(i, x0, x1,p))
        x0 = x1
        x1 = p
    return p