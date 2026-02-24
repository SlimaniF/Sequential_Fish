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
<<<<<<< HEAD
from math import comb

def compute_unique_position_expectancy(a1,V) :
    """
    Compute expected of unique position that will be occupied by a1 single molecule in V positions.
    """
    
    a1_unique = V*(1-(1-(1/V))**a1)
=======
from typing import List, Union
from pandas import Series

def p(abundancy : int, volume: int) :
    """
    The probability that a voxel m is occupied with at least one molecule.
    """
    if isinstance(volume, int) :
        if volume == 0 : return np.NaN #No positions available.

    return 1 - np.power(1-1/volume,abundancy)

def q(abundancy : int, volume : int) :
    """
    The probability that voxel m and l != m are both occupied with at least one molecule i
    """
    if isinstance(volume, int) :
        if volume == 0 : return np.NaN #No positions available.

    return 1-2*np.power((1-1/volume),abundancy) + np.power((1-2/volume),abundancy)

def c(abundancy : int, volume : int) :
    """
    The covariance of the occupancy of two different voxel by the same distribution i.
    """

    return q(abundancy, volume) - np.power(p(abundancy, volume), 2)


def compute_unique_position_expectancy(a1,V) :
    """
    Compute expectancy of unique position that will be occupied by at least one molecule with a distribution of a1 single molecule randomly spread in V positions.
    """
    if V == 0 : return np.NaN #No positions available.
    
    a1_unique = V*p(a1,V)
>>>>>>> Modification_drift_correction
    return a1_unique

def compute_self_colocalization_expectancy(a1, V) :
    """
    Compute expected number of self-colocalization which correspond to the expected number of detection that didn't discover a new pixel amongst the V pixels.
    """
<<<<<<< HEAD
    
    p = 1-(1-1/V)**a1 #Probability that a new molecule occupies an already taken position
    expected_number_selfcolocalisation = a1- V*p # total number of detection - expected number of detection landing in an already taken position = number of position that landed in a new pixel.
    
    return expected_number_selfcolocalisation

def compute_colocalization_count_expectancy(a1_unique, a2, V) :
=======

    return a1 * p(a1 -1,V)


def compute_self_colocalization_std(a1, V) :
    """
    Compute standard deviation number of self-colocalization which correspond to the expected number of detection that didn't discover a new pixel amongst the V pixels.
    """
    if V == 0 : return np.NaN #No positions available.
    p1 = p(a1-1,V)
    expected_number_selfcolocalisation = np.sqrt(a1*p1*(1-p1))
    
    return expected_number_selfcolocalisation

def compute_colocalization_probability(a1,V) :
    """
    Statistically colocalization probabilty is colocalization rate.  
    This is compute as number of **occupied position** / **available positions**.  
    
    """

    
    return p(a1,V)

def compute_colocalization_count_expectancy(a1, a2, V) :
>>>>>>> Modification_drift_correction
    """
    Compute the expected number of colocalization events which corresponds to a binomial law of sucess probability of a2 picking a position occupied by the a1_unique particules amongst V positions.
    Expectancy = np
    
    """
    
<<<<<<< HEAD
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
=======
    
    coloc_count = a2 * compute_colocalization_probability(a1,V)
    
    return coloc_count


def compute_colocalization_count_std(a1, a2, V) :
>>>>>>> Modification_drift_correction
    """
    Compute standard deviation of expected number of colocalization events which corresponds to a binomial law of sucess probability of a2 picking a position occupied by the a1_unique particules amongst V positions.
    Var = np(1-p)
    """
    
<<<<<<< HEAD
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
=======
    p2 = p(a2,V)
    c2 = c(a2,V)

    variance = a1*((1-1/V) * (p2*(1-p2) - c2) + c2)

    return np.sqrt(variance)

def Ncolocalization_expectancy(abundancies:List[Union[int,Series]], volume:Union[int,Series]) :
    """
    Expectancy deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
    """
    if isinstance(volume, int) and isinstance(abundancies,list):
        if all([isinstance(a,int) for a in abundancies]) :
            p_combination = np.prod([p(abundancy, volume) for abundancy in abundancies])
        else :
            raise TypeError("Volume is int but non int values passed in abundancies")
    
    elif isinstance(volume, Series) and isinstance(abundancies,list) :
        if all([isinstance(a,Series) for a in abundancies]) :
            p_combination = np.prod([p(abundancy, volume) for abundancy in abundancies], axis=0)
        else :
            raise TypeError("Volume is pandas Series but abundancies is not a list of Series.")

    else :
       raise TypeError("Either pass abundancies as List[int] and volume as int or abudancies as List[Series] and volume as Series") 

    return (volume*p_combination).rename("expectancy")

def Ncolocalization_std(abundancies:List[Union[int,Series]], volume:Union[int,Series]):
    """
    Standard deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
    """
    if isinstance(volume, int) and isinstance(abundancies,list):
        if all([isinstance(a,int) for a in abundancies]) :
            p_combination = np.prod([p(abundancy, volume) for abundancy in abundancies])
            q_combination = np.prod([q(abundancy, volume) for abundancy in abundancies])
        else :
            raise TypeError("Volume is int but non int values passed in abundancies")
    elif isinstance(volume, Series) and isinstance(abundancies,list) :
        if all([isinstance(a,Series) for a in abundancies]) :
            p_combination = np.prod([p(abundancy, volume) for abundancy in abundancies], axis=0)
            q_combination = np.prod([q(abundancy, volume) for abundancy in abundancies], axis=0)
        else :
            raise TypeError("Volume is pandas Series but abundancies is not a list of Series.")



    variance = volume*p_combination*(1-p_combination) + volume*(volume-1)*(q_combination - np.power(p_combination,2))
    
    return np.sqrt(variance).rename("std")
>>>>>>> Modification_drift_correction
