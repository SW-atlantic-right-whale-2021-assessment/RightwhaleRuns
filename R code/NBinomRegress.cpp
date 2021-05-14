// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Read data
    DATA_MATRIX(X);
    DATA_MATRIX(Xd);
    DATA_VECTOR(yobs);

    // Parameters
    PARAMETER_VECTOR(beta); // Expected value parameters
    PARAMETER_VECTOR(betad); // dispersion parameters

    // Variables
    vector<Type> eta = X*beta; // Expected value parameters
    vector<Type> etad = Xd*betad; // Dispersion

    // Run estimation
    for (int i=0; i < yobs.size(); i++){
        s1 = eta(i);          // log(mu)
        s2 = 2. * s1 - etad(i) ;  // log(var - mu)
        jnll -= dnbinom_robust(yobs(i), s1, s2, true);
    }
    return jnll;
}