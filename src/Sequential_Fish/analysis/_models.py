"""
Submodule containing implementation of statistical/probabilistic models
"""

"""
1. Colocalisation model --> multicolocalization_statistics.ipynb

    We want to modelise the colocalization dynamics in case of random positioning of single molecules in the availables pixel.
    In this notation we try the colocalization between target (rna) 1 with target 2, meaning we want to know how many single molecule from
    target 1 distribution are expected to co-localize with at least one single molecule of target 2.

    * V : the total number of available positions (i.e pixel) representing the volume of a cell.
    * a1 : abundancy of target 1, i.e number of unique molecule of target 1.
    * a2 : abundancy of target 2, i.e number of unique molecule of target 2.
    * ai_unique : number of unique positions occupied by at least one single molecule of target i.

"""

import numpy as np
from math import comb

def compute_unique_position_expectancy(a1,V) :
    """
    Compute expected of unique position that will be occupied by a1 single molecule in V positions.
    """
    
    a1_unique = V*(1-(1-(1/V))**a1)
    return a1_unique

def compute_self_colocalization_expectancy(a1, V) :
    """
    Compute expected number of self-colocalization which correspond to the expected number of detection that didn't discover a new pixel amongst the V pixels.
    """
    
    p = 1-(1-1/V)**a1 #Probability that a new molecule occupies an already taken position
    expected_number_selfcolocalisation = a1- V*p # total number of detection - expected number of detection landing in an already taken position = number of position that landed in a new pixel.
    
    return expected_number_selfcolocalisation

def compute_colocalization_count_expectancy(a1_unique, a2, V) :
    """
    Compute the expected number of colocalization events which corresponds to a binomial law of sucess probability of a2 picking a position occupied by the a1_unique particules amongst V positions.
    Expectancy = np
    
    """
    
    if V == 0 : return np.nan #No positions available.
    
    coloc_count = a2 * a1_unique/V
    
    return coloc_count

def compute_colocalization_probability(a1_unique,V) :
    """
    Statistically colocalization probabilty is colocalization rate.  
    This is compute as number of **occupied position** / **available positions**.  
    
    """
    
    p = a1_unique / V # probablity to pick a red ball if there are a1_unique red balls amongst V balls.
    
    return p
    

def compute_colocalization_count_std(a1_unique, a2, V) :
    """
    Compute standard deviation of expected number of colocalization events which corresponds to a binomial law of sucess probability of a2 picking a position occupied by the a1_unique particules amongst V positions.
    Var = np(1-p)
    """
    
    if V == 0 : return np.nan #No positions available.
    
    std = np.sqrt(a2 * a1_unique/V * (1-a1_unique/V))
    
    return std

def compute_unique_pair_expectancy(a1_unique, a2, V) :
    """
    Compute expected number of unique pair formed between a2 singles and a1_unique singles when the number of positions available is V.
    """
    if V == 0 : return np.nan
    
    unique_pair_count = a1_unique*(1-(1-1/V)**a2)
    
    return unique_pair_count


def compute_unique_pair_std(a1_unique, a2, V) :
    """
    Compute expected standard deviation of number of unique pair formed between a2 singles and a1_unique singles when the number of positions available is V.
    
    We use the formula of variance for a sum of random variables Ii.
    var = iSum(Var(I)) + 2*iSum(Cov(I))
    
    First term is called Exp
    second Cov
    
    Exp = Sum(p*(1-p)) (Bernoulli law)
     --> a1_unique *(p(1-p))
    
    Cov = E(Ii and Ij) - E(Ii)*E(Ij) (called Eij, Ei and Ej)
    
    """
    if V == 0 : return np.nan
    
    p = 1 - (1-1/V)**a2 # probability to draw a specific location at least once
    Exp = a1_unique*p*(1-p) # Variance of bernouli law
    Cov = (1-2/V)**a2 - (1-1/V)**(2*a2) #Simplified expression for cov
    
    var = Exp + 2*comb(round(a1_unique),2)*Cov
    
    return np.sqrt(var)